import sys

INPUT_FILE = "in.txt"
OUTPUT_FILE = "out.txt"

def get_idx(num, base=2, exponent=4)->list[int]:
    idx = []
    for i in range(exponent):
        idx.append(num % base)
        num //= base
    return idx

class CFD:
    def __init__(self, parents, row_cfd):
        self.parents = parents
        self.cfd = row_cfd

    def query(self, result, evidence):
        cfd_idx = sum([2**(len(self.parents)-i-1) * evidence[parent] for i, parent in enumerate(self.parents)])
        return self.cfd[cfd_idx] if result == 1 else 1 - self.cfd[cfd_idx]

class RV:
    def __init__(self, x: int):
        self.x = x
        self.parents: list[RV] = []
        self.children: list[RV] = []

    def set_cfd(self, _cfd: CFD):
        self.cfd = _cfd

    def add_child(self, child):
        self.children.append(child)

    def set_parents(self, parents):
        self.parents = parents
    
    def __hash__(self):
        return hash(self.x)
    def __equal__(self, rhs):
        return self.x == rhs.val()
    def __lt__(self, other):
        return self.x < other.val()
    def val(self):
        return self.x

class BN:
    def __init__(self, n_rvs):
        self.X: list[RV] = [RV(x) for x in range(n_rvs)]
    
    # add parents and cfd
    def add_relation(self, row_node, row_parents, row_cfd):
        node = self.X[row_node]
        parents = [self.X[parent] for parent in row_parents]
        cfd = CFD(parents, row_cfd)

        node.set_parents(parents)
        node.set_cfd(cfd)
        for parent in parents:
            parent.add_child(node)
    
    def MarkovBlanket(self, node: RV, exclude=[]):
        mb:set[RV] = set()
        for child in node.children:
            if child in exclude: continue
            mb.add(child)
            for parent in child.parents:
                mb.add(parent)
        mb = mb.difference(exclude)
        if node not in exclude:
            for parent in node.parents: mb.add(parent)
        return mb.difference([node])

    def query(self, rv: int, E: dict):
        exclude = set()
        eliminated = set()
        message_val = set()
        message = dict()

        # for debugging print message
        def print_message():
            for k, v in message.items():
                print (*[f"{n.val()}: {p}" for n, p in k], "->", v)

        # for looping every node and eliminate (except rv)
        for x in self.X:
            print(x.val())
            if x.val() == rv: continue

            children = set(x.children).difference(exclude)
            MB_x = self.MarkovBlanket(x, exclude)
            if x in message_val: MB_x = MB_x.union(message_val).difference([x])
            new_message = dict()

            print("chld: ", *[ii.val() for ii in children])
            print("MB_X: ", *[ii.val() for ii in MB_x])
            
            # for looping every possible values of MB(x)
            for i in range(2**len(MB_x)):
                # val = {element of MB_x : value (0 or 1)}
                val = {_mb: _v for _mb, _v in zip(MB_x, get_idx(i, exponent=len(MB_x)))}
                print("val: ", *[f"{_mb.val()}: {_v}" for _mb, _v in val.items()])
                
                # check if evidence is setted correctly
                invalid=False
                for key in E.keys():
                    if key in val and val[key] != E[key]: invalid=True
                if invalid: continue
                
                # calculate probability
                prob = [1, 1]
                if x not in exclude:
                    prob = [x.cfd.query(0, val), x.cfd.query(1, val)]
                for child in children:
                    print(prob)
                    for x_val in range(2):
                        val[x] = x_val
                        prob[x_val] *= child.cfd.query(val[child], val)
                print(prob)
                if len(children) > 0: del val[x]

                if x in message_val:
                    for x_val in range(2):
                        if x in E and E[x] != x_val: continue
                        query = {m:val[m] for m in message_val if m!=x}
                        query[x] = x_val
                        prob[x_val] *= message[tuple(sorted(query.items()))]
                new_message[tuple(sorted(val.items()))] = prob[E[x]] if x in E else sum(prob)
            
            message = new_message
            print_message()
            message_val = {mv for mv,_ in list(message.keys())[0]}
            eliminated.add(x)
            exclude.add(x)
            exclude = exclude.union(x.children)
        
        print(message)
            

def get_input():
    N = int(input())
    bn = BN(N)

    for i in range(N):
        parent_info = list(map(lambda x: int(x)-1, input().split(" ")))[1:]
        cfd_info = list(map(float, input().split(" ")))
        bn.add_relation(i, parent_info, cfd_info)
    
    return bn
    



if __name__ == "__main__":
    input_file = open(INPUT_FILE, 'r')
    output_file = open(OUTPUT_FILE, 'w')

    sys.stdin = input_file
    # sys.stdout = output_file

    # get input and construct BN
    bn = get_input()
    for x in bn.X:
        print(x.val(), x.cfd.cfd)
        print("parents: ", [i.val() for i in x.parents])
        print("children: ", [i.val() for i in x.children])

    bn.query(2, {bn.X[0]:1,bn.X[1]:0})

    input_file.close()
    output_file.close()


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

class Message:
    def __init__(self, associated_node):
        self.associated_node = associated_node
        self.prob = dict()
    
    def add_prob(self, val, prob):
        self.prob[tuple(sorted(val.items()))] = prob
    def query(self, val):
        _val = {_k:_v for _k, _v in val.items() if _k in self.associated_node}
        return self.prob[tuple(sorted(_val.items()))]
    def __contains__(self, x):
        return x in self.associated_node
    

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

    def variable_elimination(self, rv: list[int], E: dict):
        exclude = set()
        all_messages: list[Message] = []
        
        # for _rv in rv:
        #     if len(self.X[_rv].parents) == 0:
        #         new_message = Message({self.X[_rv]})
        #         new_message.add_prob({self.X[_rv]: 0}, self.X[_rv].cfd.query(0, dict()))
        #         new_message.add_prob({self.X[_rv]: 1}, self.X[_rv].cfd.query(1, dict()))
        #         all_messages.append(new_message)

        # all_messages printer for debugging
        def _print_messages():
            for message in all_messages:
                for k, v in message.prob.items():
                    print (*[f"{n.val()}: {p}" for n, p in k], "->", v, file=sys.stderr)

        # for looping every node in BN except "rv" node (to eliminate them)
        for x in self.X:
            print(x.val(), file=sys.stderr)
            _print_messages()
            
            ## candidates
            # candidates are random variables that its cfd should be used
            # we remove variables in exclude because the information about them are in all_messages
            ## associated_node
            # associated_node variable contains every node which needs to construct new message
            candidates = set(x.children).difference(exclude)
            associated_node = set([pa for can in candidates for pa in can.parents] + [can for can in candidates])
            if x not in exclude: associated_node = associated_node.union(set(x.parents))
            
            # pick messages that related with rv x and remove them from all_messages
            # Make another new message(new_message) that carries the probability which will be calculated this iteration
            messages: list[Message] = []
            _messages = all_messages.copy()
            for m in _messages:
                if x in m:
                    messages.append(m)
                    all_messages.remove(m)
            for message in messages:
                associated_node = associated_node.union(message.associated_node)
            associated_node = associated_node.difference([x])
            new_message = Message(associated_node if x.val() not in rv else associated_node.union([x]))
            
            # for looping to assign every possible value to associated_node
            for _val in range(2**len(associated_node)):
                # val is possible values the associated_node can have
                # {element of associated_node: value (0 or 1)}
                val = {_node: _v for _node, _v in zip(associated_node, get_idx(_val, exponent=len(associated_node)))}


                # if associated_node includes some evidence nodes, check if the val is setted correctly.
                # if not, continue the loop
                invalid = False
                for key in E.keys():
                    if key in val and val[key] != E[key]: invalid=True
                if invalid: continue

                # calculate probability for given "val"
                prob = [1, 1]
                if x not in exclude:
                    prob = [x.cfd.query(0, val), x.cfd.query(1, val)]

                for x_val in range(2):
                    if x in E and E[x] != x_val: continue
                    # initialize val and prob
                    val[x] = x_val

                    # multipies probability from cfd and messages
                    for candidate in candidates:
                        prob[x_val] *= candidate.cfd.query(val[candidate], val)
                    for m in messages:
                        prob[x_val] *= m.query(val)

                # if x is from rv, save probability separately
                if x.val() in rv:
                    val[x]=0
                    new_message.add_prob(val, prob[0])
                    val[x]=1
                    new_message.add_prob(val, prob[1])
                else:
                    del val[x]
                    final_prob = prob[E[x]] if x in E else sum(prob)
                    new_message.add_prob(val, final_prob)
            
            # update messages and add x and its children to exclude
            # including them means that their CFD wouldn't be used
            all_messages.append(new_message)
            exclude.add(x)
            exclude = exclude.union(x.children)
        _print_messages()
        return all_messages
    
    def query(self, _X, _E):
        E = {self.X[_k]: _v for _k, _v in _E.items()}
        P_X_E = self.variable_elimination(_X, E)
        P_E = self.variable_elimination(_E.keys(), dict())
        return P_X_E[0].query({self.X[_x]: 1 for _x in _X}) / P_E[0].query(E)

            

def get_input():
    N = int(input())
    bn = BN(N)

    # construct BN
    for i in range(N):
        parent_info = list(map(lambda x: int(x)-1, input().split(" ")))[1:]
        cfd_info = list(map(float, input().split(" ")))
        bn.add_relation(i, parent_info, cfd_info)
    
    # get input about X and evidence
    X, n = tuple(map(int, input().split(" ")))
    X-=1
    E = dict()
    for i in range(n):
        node, val = tuple(map(int, input().split(" ")))
        E[node-1] = val
    
    return bn, X, E
    



if __name__ == "__main__":
    input_file = open(INPUT_FILE, 'r')
    output_file = open(OUTPUT_FILE, 'w')

    sys.stdin = input_file
    sys.stdout = output_file

    # get input and construct BN
    bn, X, E = get_input()
    print(bn.query([X], E))

    input_file.close()
    output_file.close()


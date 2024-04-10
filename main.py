import sys

INPUT_FILE = "in.txt"
OUTPUT_FILE = "out.txt"

class CFD:
    def __init__(self, parents, row_cfd):
        self.parents = parents
        self.cfd = row_cfd

class RV:
    def __init__(self, x: int, n_space:int = 2):
        self.x = x
        self.n_space = n_space
        self.parents: list[RV] = []
        self.children: list[RV] = []

    def set_cfd(self, _cfd):
        self.cfd = _cfd

    def add_child(self, child):
        self.children.append(child)

    def set_parents(self, parents):
        self.parents = parents
    
    def __hash__(self):
        return hash(self.x)
    def __equal__(self, rhs):
        return self.x == rhs.val()
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

def get_input():
    N = int(input())
    bn = BN(N)

    for i in range(N):
        parent_info = list(map(lambda x: int(x)-1, input().split(" ")))[1:]
        cfd_info = list(map(int, input().split(" ")))
        bn.add_relation(i, parent_info, cfd_info)
    
    return bn
    



if __name__ == "__main__":
    input_file = open(INPUT_FILE, 'r')
    output_file = open(OUTPUT_FILE, 'w')

    sys.stdin = input_file
    sys.stdout = output_file

    # get input and construct BN
    bn = get_input()



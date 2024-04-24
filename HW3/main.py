import sys 

LM = 1000

N = 0
alpha, beta, gamma = 0., 0., 0.
graph = []
neighbors = [0 for _ in range(LM)]
dfs, visited = [], [0 for _ in range(LM)]
message = [[() for _ in range(LM)] for _ in range(LM)]

class Node:
    def __init__(self, n):
        self.n = n
        self.neighbors = []
    
    def add_neighbor(self, neighbor):
        i = 0
        while i < len(self.neighbors):
            if self.neighbors[i].n > neighbor.n: break
            i+=1
        self.neighbors.insert(i, neighbor)
    
    def __eq__(self, node):
        return node.n == self.n


def factor(i, j, X_i, X_j):
    if i ==j:
        return 2 * (alpha * i * X_i + gamma)
    elif i < j:
        return alpha * i * X_i + beta * j * X_j + gamma
    else:
        return alpha * j * X_j + beta * i * X_i + gamma



def get_input():
    global N, alpha, beta, gamma, graph

    N = int(input())
    graph = [Node(i) for i in range(N)]
    for _ in range(N-1):
        a, b = tuple(map(int, input().split(" ")))
        neighbors[a] += 1
        neighbors[b] += 1
        graph[a].add_neighbor(graph[b])
        graph[b].add_neighbor(graph[a])
    alpha, beta, gamma = tuple(map(float, input().split(" ")))


# for debugging (printing graph structure)
def print_graph():
    for i in range(N):
        print(i, end="   ")
    print()
    
    for i in range(N-1):
        o_stream = [" " for _ in range(4 * N - 3)]
        for node in graph:
            pos = node.n * 4
            for neighbor in node.neighbors:
                diff = neighbor.n - node.n
                if diff <= i and diff >= -i: continue

                sign = pos + (i+1 if diff > 0 else -i-1)
                if (o_stream[sign] == "/" and diff > 0) or (o_stream[sign] == "\\" and diff < 0): o_stream[sign] = "X"
                else: o_stream[sign] = "\\" if diff > 0 else "/"

                if diff == i+1:
                    for j in range(sign+1, sign+2*(i+1)): o_stream[j] = "_"
        print(''.join(o_stream))


def run_forward_dfs():
    while len(dfs) != 0:
        cur_node = dfs.pop(0)

        # product of input message
        product = [1, 1]
        # passing message to this value
        dest = None

        # for looping neighbors and calculate message from neighbors
        for neighbor in cur_node.neighbors:
            # if neighbor gave message to cur_node, calculate message
            if len(message[neighbor.n][cur_node.n]) != 0:
                product[0] *= message[neighbor.n][cur_node.n][0]
                product[1] *= message[neighbor.n][cur_node.n][1]
            else: dest = neighbor
        
        # calculate message from cur_node to dest
        message_ji = [0, 0]
        i, j = dest.n, cur_node.n
        for _i in range(2):
            X_i = 2*_i-1
            for _j in range(2):
                X_j = 2*_j-1
                message_ji[_i] += factor(j, j, X_j, X_j) * factor(i, j, X_i, X_j) * product[_j]
        message[j][i] = tuple(message_ji)

        # if dest is visited one less than its neighbors, it is prepared to calculate message
        visited[dest.n] += 1
        if visited[dest.n] == neighbors[dest.n]-1: dfs.append(dest)

def run_backward_dfs():
    while len(dfs) != 0:
        cur_node = dfs.pop(0)

        # for looping neighbors and calculate message from neighbors
        for neighbor in cur_node.neighbors:
            if visited[neighbor.n]: continue
            # product of input message
            product = [1, 1]
            # passing message to this value
            dest = neighbor

            # calculate products
            for _neighbor in cur_node.neighbors:
                if _neighbor == dest: continue
                product[0] *= message[_neighbor.n][cur_node.n][0]
                product[1] *= message[_neighbor.n][cur_node.n][0]

            # calculate message from cur_node to dest
            message_ji = [0, 0]
            i, j = dest.n, cur_node.n
            for _i in range(2):
                X_i = 2*_i-1
                for _j in range(2):
                    X_j = 2*_j-1
                    message_ji[_i] += factor(j, j, X_j, X_j) * factor(i, j, X_i, X_j) * product[_j]
            message[j][i] = tuple(message_ji)

            # check visited after pushing it to dfs queue
            visited[dest.n] = True
            dfs.append(dest)


def message_passing():
    global visited
    # forward
    # initialize dfs queue
    root = Node(-1)
    for i in range(N):
        if neighbors[i] == 1:
            if root == Node(-1): root = graph[i]
            else: dfs.append(graph[i])
    run_forward_dfs()
    
    # backward
    # initialize dfs queue
    dfs.append(root)
    visited = [False for _ in range(N)]
    visited[0] = True
    run_backward_dfs()
    

if __name__ == "__main__":
    input_file = open("input.txt", 'r')
    output_file = open("output.txt", 'w')

    sys.stdin = input_file
    # sys.stdout = output_file

    get_input()

    print("N: ", N)
    print("a, b, c: ", alpha, beta, gamma)
    print_graph()

    message_passing()

    for val in range(2):
        for node in graph:
            for neighbor in node.neighbors:
                if node.n < neighbor.n:
                    print(message[node.n][neighbor.n][val], end=" ")
        print()
        for node in graph:
            for neighbor in node.neighbors:
                if node.n < neighbor.n:
                    print(message[neighbor.n][node.n][val], end=" ")
        print()


    input_file.close()
    output_file.close()
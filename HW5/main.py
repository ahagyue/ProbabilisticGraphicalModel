import numpy as np
import math
from scipy.special import polygamma
import matplotlib.pyplot as plt

class Documents:
    def __init__(self):
        self.M = 0
        self.document = dict()
        self.Nd = dict()
        self.L = dict()
    
    def add_word(self, _doc_id, w, frequence):
        doc_id = _doc_id - 1
        if doc_id in self.document:
            self.document[doc_id].append((w, frequence))
            self.Nd[doc_id] += 1
            self.L[doc_id] += frequence
        else:
            self.M += 1
            self.document[doc_id] = [(w, frequence)]
            self.Nd[doc_id] = 1
            self.L[doc_id] = frequence
    
    def print_heads(self):
        for d in range(10):
            for j, w in enumerate(self.document[d]):
                print(d, w)
                if j == 3: break
            print("...")

def digamma(X):
    op1 = 0.00416667 * 1 / ((X + 6)**2) - 0.003968
    op2 = op1 * 1 / ((X + 6) ** 2) + 0.0083333
    op3 = op2 * 1 / ((X + 6) ** 2) - 0.0833333
    op4 = op3 * 1 / ((X + 6) ** 2) + np.log(X) - 1/(2*X)
    return op4 - 1 / (X - 1) - 1 / (X - 2) - 1 / (X - 3) - 1 / (X - 4) - 1 / (X - 5) - 1 / (X - 6)

def calc_likelihood(document, N, K, alpha, beta, phi, gamma, digam):
    gam_sum = np.sum(gamma)
    digam_sum = np.sum(digam)
    like = math.lgamma(alpha * K) - K * math.lgamma(alpha) - math.lgamma(gam_sum)
    
    for k in range(K):
        like += (alpha - 1) * (digam[k] - digam_sum) + math.lgamma(gamma[k]) - (gamma[k] - 1) * (digam[k] - digam_sum)
        for n in range(N):
            if phi[n, k] > 0:
                like += document[n][1] * (phi[n, k] * (digam[k] - digam_sum - math.log(phi[n][k]) + beta[k, document[n][0] - 1]))
    return like


def LDA(K, documents):
    def logsum(a, b):
        M = a if a > b else b
        m = a if a <= b else b
        return M + math.log(1 + math.exp(m - M))

    V = 61188
    M = documents.M

    alpha   = 50 / K
    beta    = np.zeros((K, V))   # K X V
    
    like = 0
    likelihood_prev = 0
    like_hist = []
    converged = False

    MAX_IT = 20
    EM_MAX_IT = 50
    DIFF_LIM = 1e-6
    EM_LIM = 1e-3

    iterations = 0
    while not converged:
        iterations += 1
        print(f"{iterations}th EM iteration")

        classword = np.zeros((K, V))
        classtotal = np.zeros(K)
        ndocs = 0
        alpha_ss = 0
        likelihood = 0

        # E-step
        for d in range(100):
            if d % 200 == 0: print(f"{d}th E step document")

            oldlike = 0

            N = documents.Nd[d]
            L = documents.L[d]
            rand    = np.random.rand(N, K)
            phi     = rand / rand.sum(axis=0)                       # N X K
            gamma   = np.ones(K) * alpha + np.sum(phi, axis=0)      # K
            digam   = polygamma(1, gamma)                           # K
            oldphi  = np.zeros(K)                                   # K

            E_it = 0
            while True:
                E_it += 1
                for n in range(N):
                    phisum = 0
                    for k in range(K):
                        oldphi[k] = phi[n][k]
                        phi[n][k] = digam[k] + beta[k, documents.document[d][n][0] - 1]
                        phisum = phi[n][k] if k == 0 else logsum(phisum, phi[n][k])
                    
                    for k in range(K):
                        phi[n, k] = math.exp(phi[n, k] - phisum)
                        gamma[k] = gamma[k] + documents.document[d][n][1] * (phi[n, k] - oldphi[k])
                        digam[k] = polygamma(1, gamma[k])

                like = calc_likelihood(documents.document[d], N, K, alpha, beta, phi, gamma, digam)
                like_diff = (oldlike - like) / oldlike
                oldlike = like

                if (like_diff < DIFF_LIM or E_it > MAX_IT):
                    break
            
            gam_sum = np.sum(gamma)
            alpha_ss += sum(polygamma(1, gamma)) - K * polygamma(1, gam_sum)

            for n in range(N):
                for k in range(K):
                    classword[k, documents.document[d][n][0] - 1] += documents.document[d][n][1] * phi[n, k]
                    classtotal[k] += documents.document[d][n][1] * phi[n, k]

            ndocs += 1
            likelihood += like
            if like == np.nan: 
                print(d, " occured nan")
                break
        
        # M-step
        print("M step")
        for k in range(K):
            for w in range(V):
                if classword[k, w] > 0:
                    beta[k, w] = math.log(classword[k, w]) - math.log(classtotal[k])
                else:
                    beta[k, w] = -100
        
        # alpha estimation
        alpha_iter = 0
        a_init = 100
        a_log = math.log(a_init)
        while True:
            alpha_iter += 1
            a = math.exp(a_log)

            df = ndocs * (K * polygamma(1, K * a) - K * polygamma(1, a)) + alpha_ss
            d2f = ndocs * (K * K * polygamma(2, (K * a)) - K * polygamma(2, a))
            a_log -= df / (d2f * a + df)
            print("alpha update: ", math.exp(a_log))
            if (math.fabs(df) < 1e-5 or alpha_iter > 100):
                break
        alpha = math.exp(a_log)

        likelihood_diff = (likelihood_prev - likelihood) / likelihood_prev
        likelihood_prev = likelihood
        like_hist.append(likelihood)

        print(likelihood_diff, likelihood)
        if likelihood_diff < 0: MAX_IT *=2
        if ((likelihood_diff < EM_LIM and likelihood_diff > 0) or iterations > EM_MAX_IT) and iterations > 2:
            break
    
    plt.plot(np.arange(len(like_hist)), like_hist)
    plt.show()

    return alpha, beta


def EDA(prob, K, n_most = 16):
    vocabs = []
    with open("./data/20newsgroup/vocabulary.txt") as f:
        while True:
            vocab = f.readline()[:-1]
            if not vocab: break
            vocabs.append(vocab)
    
    ind = prob.argsort()[:, -n_most:][:, ::-1]

    for k in range(K):
        print("%-30s" % f"{k+1} topc", end="")
    print()

    for n in range(n_most):
        for k in range(K):
            print("%-30s" % f"{vocabs[ind[k, n]]}: {round(math.exp(prob[k, ind[k, n]]), 2)}", end="")
        print()


if __name__ == "__main__":
    K = 10

    train_documents = Documents()
    test_documents = Documents()

    with open("./data/20newsgroup/train.data") as f:
        while True:
            line = f.readline()
            if not line: break
            train_documents.add_word(*map(int, line.split()))
    
    with open("./data/20newsgroup/test.data") as f:
        while True:
            line = f.readline()
            if not line: break
            test_documents.add_word(*map(int, line.split()))

    alpha, beta = LDA(K, train_documents)

    EDA(beta, K, 1000)
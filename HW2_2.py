probability_set = {
    (0, 0, 0, 0): 1/8,
    (1, 1, 0, 0): 1/8,
    (1, 1, 1, 0): 1/4,
    (0, 1, 0, 1): 1/4,
    (1, 0, 1, 1): 1/4,
}

def get_idx(num, base=2, exponent=4):
    idx = []
    for i in range(exponent):
        idx.append(num % base)
        num //= base
    return idx

# marginalization of all cases
for every_cases in range(3**4):
    idx = get_idx(every_cases, base=3)

    prob = 0
    marginals = [i for i in range(4) if idx[i] == 2]
    for marginal_rv in range(2**len(marginals)):
        rv_vals = get_idx(marginal_rv, exponent=len(marginals))
        joint_rv = [0 if i in marginals else idx[i] for i in range(4)]
        for rv_val_idx, marginal_idx in enumerate(marginals): joint_rv[marginal_idx] = rv_vals[rv_val_idx]
        joint_rv = tuple(joint_rv)
        
        prob += 0 if not joint_rv in probability_set else probability_set[joint_rv]
    probability_set[tuple(idx)] = prob


independence = []

for A in range(3):
    for B in range(A+1, 4):
        _E = [i for i in range(4) if i!=A and i!=B]
        E = [[_E[j] for j, k in enumerate(get_idx(i, exponent=2)) if k==1] for i in range(2**2)]
        for evidence in E:
            is_independent = True
            ev_cases = []
            for i in range(2**len(evidence)):
                case = [2, 2, 2, 2]
                for j, k in zip(evidence, get_idx(i, exponent=len(evidence))):
                    case[j] = k
                ev_cases.append(case)
            for ev_case in ev_cases:
                for i in range(2):
                    for j in range(2):
                        AE = [i for i in ev_case]
                        BE = [i for i in ev_case]
                        ABE = [i for i in ev_case]
                        AE[A] = i
                        BE[B] = j
                        ABE[A] = i
                        ABE[B] = j
                        P_E = probability_set[tuple(ev_case)]
                        if P_E == 0: continue
                        P_AE = probability_set[tuple(AE)]
                        P_BE = probability_set[tuple(BE)]
                        P_ABE = probability_set[tuple(ABE)]

                        if (P_ABE / P_E) != (P_AE*P_BE/(P_E**2)):
                            is_independent = False
                        
            if is_independent:
                independence.append(((A, B), evidence))
print(independence)
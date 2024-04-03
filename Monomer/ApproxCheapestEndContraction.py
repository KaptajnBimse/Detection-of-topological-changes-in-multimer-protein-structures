import numpy as np
from PriceEstEndContraction import PriceEstEndContraction

def ApproxCheapestEndContraction(M, len):
    M = M[np.argsort(M[:, 4])]
    Nbr = M.shape[0]
    B = np.zeros((Nbr+1, 3))
    B[0, 0] = 1
    B[0, 1] = min(np.concatenate((M[:, 4], [len])))
    B[1:, 0] = M[:, 5]
    for i in range(Nbr):
        B[i+1, 1] = min(np.concatenate((M[i+1:, 4], [len])))
    B[:, 2] = PriceEstEndContraction(B[:, 0]-1) + PriceEstEndContraction(len-B[:, 1])
    I = np.argmin(B[:, 2])
    b = B[I, 1]
    a = B[I, 0]
    return a, b
import numpy as np
from intersection_origo_triangle_line_segment import intersection_origo_triangle_line_segment
import d_points2line as dpl

def IsContractableType1ReparametrizationParallel(M, M0, M1, i, P, P1, maxlen):
    FindNumberOfOmega1_2Obstructions = 0
    i = int(i)-1 # i is 1-indexed in the original code 
    sav = M0[i, 7]
    P = ((1 - sav) * P + sav * P1).T
    mint1 = M0[i, 4]
    maxt1 = M0[i, 3]
    leng1 = maxt1 - mint1
    mint2 = M1[i, 4]
    maxt2 = M1[i, 3]
    leng2 = maxt2 - mint2
    looplength = max(leng1, leng2)
    if looplength > maxlen:
        return [0, 0]

    mint = M[i, 4]
    maxt = M[i, 3]
    avt = (mint + maxt) / 2
    n1av = int(np.floor(avt))
    tav = avt - n1av

    n1 = int(np.floor(mint))
    n2 = int(np.ceil(maxt))
    a = mint - n1
    b = maxt - np.floor(M[i, 3])

    pts = np.column_stack(((1 - a) * P[:, n1-1] + a * P[:, n1], P[:, n1 + 0:(n2-1)], (1 - b) * P[:, n2-2] + b * P[:, n2-1]))
    if np.sum((pts[:, 0] - pts[:, -1]) ** 2) > 10**(-15):
        pointdistance = np.sum((pts[:, 0] - pts[:, -1]) ** 2) ** 0.5
        print('WARNINGNoIntersection distance', pointdistance)
        return [0, 0]

    center = np.sum(pts[:, 0:-1], axis=1) / (pts.shape[1] - 1)
    pts = pts - np.tile(center.reshape(-1, 1), (1, pts.shape[1]))
    rdisk = np.max(np.sum(pts ** 2, axis=0) ** 0.5)
    pmidt = (1 - tav) * P[:, n1av] + tav * P[:, n1av + 1] - center

    NbrTriangles = pts.shape[1] - 1
    P = P - np.tile(center.reshape(-1, 1), (1, P.shape[1]))
    Lstart = P[:, np.r_[0:n1 - 1, n2:P.shape[1] - 1]]
    Lend = P[:, np.r_[1:n1, n2 + 1:P.shape[1]]]
    Lmidt = np.sum(((Lstart + Lend) / 2) ** 2, axis=0) ** 0.5

    LineSegmentLength = np.sum((Lstart - Lend) ** 2, axis=0) ** 0.5
    ex = np.where(Lmidt <= rdisk + LineSegmentLength / 2)[0]

    Lstart = Lstart[:, ex]
    Lend = Lend[:, ex]
    NbrL = Lstart.shape[1]
    NbrIntc = 0

    if NbrL > 0:
        if FindNumberOfOmega1_2Obstructions:
            for j in range(NbrL):
                for k in range(NbrTriangles):
                    NbrIntc += intersection_origo_triangle_line_segment(pts[:, [k, k + 1]], Lstart[:, j], Lend[:, j])
        else:
            slet = np.column_stack((ex, Lmidt[ex]))
            index = np.argsort(slet[:, 1])
            for j in index:
                for k in range(NbrTriangles):
                    if intersection_origo_triangle_line_segment(pts[:, [k, k + 1]], Lstart[:, j], Lend[:, j]):
                        return [0, 0]

    if NbrIntc > 0:
        return [0, 0]

    return [np.sum(dpl.d_points2line(pts[:, 1:-1], pts[:, 0], pmidt)) * 2, looplength]
""" 
i = np.loadtxt("Test txt/IsContractableType1ReparametrizationParallel/i.txt")
M = np.loadtxt("Test txt/IsContractableType1ReparametrizationParallel/M.txt")
M0 = np.loadtxt("Test txt/IsContractableType1ReparametrizationParallel/M0.txt")
M1 = np.loadtxt("Test txt/IsContractableType1ReparametrizationParallel/M1.txt")
P = np.loadtxt("Test txt/IsContractableType1ReparametrizationParallel/P.txt")
P1 = np.loadtxt("Test txt/IsContractableType1ReparametrizationParallel/P1.txt")
maxlen = np.loadtxt("Test txt/IsContractableType1ReparametrizationParallel/maxlen.txt")
 """
#print(IsContractableType1ReparametrizationParallel(M, M0, M1, i, P, P1, maxlen))
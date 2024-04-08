import numpy as np
import intersection_origo_triangle_line_segment as iotls
import IntersectionTriangle_LineSegment as itls
import d_points2line as dp2l

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def IsContractableType2ReparametrizationParallel(M, M0, M1, i, makker, P, P1, maxlen):
    casea = "Error"
    FindNumberOfOmega1_2Obstructions = 0
    printout = 0
    printoutobstruction = 0

    leng = np.sum(np.abs(M[i, [3, 4]] - M[makker, [3, 4]]))
    lengRep = np.sum(np.maximum(np.abs(M0[i, [3, 4]] - M0[makker, [3, 4]]), np.abs(M1[i, [3, 4]] - M1[makker, [3, 4]])))

    if lengRep > maxlen:
        ud = [0, 0]
        return ud

    if printout:
        print([i, makker])
        print(M[[i, makker], :])

    ts1 = M[i, 7]
    ts2 = M[makker, 7]
    sav = (ts1 + ts2) / 2

    P0 = P
    P = ((1 - sav) * P + sav * P1).T

    Pts1 = ((1 - ts1) * P0 + ts1 * P1).T
    Pts2 = ((1 - ts2) * P0 + ts2 * P1).T
    sa = M[i, 4]
    na = np.floor(sa)
    a = sa - na
    sa2 = M[i, 3]
    na2 = np.floor(sa2)
    a2 = sa2 - na2
    sb = M[makker, 4]
    nb = np.floor(sb)
    b = sb - nb
    sb2 = M[makker, 3]
    nb2 = np.floor(sb2)
    b2 = sb2 - nb2

    La1 = (1 - a) * Pts1[:, int(na)] + a * Pts1[:, int(na+1)]
    #La1 = (1 - a) * Pts1[:, int(na)] + a * Pts1[:, int(na + 1)]
    La2 = (1 - a) * Pts2[:, int(na)] + a * Pts2[:, int(na+1)]
    La12 = (1 - a2) * Pts1[:, int(na2)] + a2 * Pts1[:, int(na2+1)]
    La22 = (1 - a2) * Pts2[:, int(na2)] + a2 * Pts2[:, int(na2+1)]
    Lb1 = (1 - b) * Pts1[:, int(nb)] + b * Pts1[:, int(nb+1)]
    Lb2 = (1 - b) * Pts2[:, int(nb)] + b * Pts2[:, int(nb+1)]
    Lb12 = (1 - b2) * Pts1[:, int(nb2)] + b2 * Pts1[:, int(nb2+1)]
    Lb22 = (1 - b2) * Pts2[:, int(nb2)] + b2 * Pts2[:, int(nb2+1)]

    mins = np.atleast_2d(np.min(M[[i, makker], 3:5], axis=0))
    maxs = np.atleast_2d(np.max(M[[i, makker], 3:5], axis=0))

    if M[i, 4] == mins[0,1]: 
        istart = i
        islut = makker
    else:
        istart = makker
        islut = i
    n1 = np.floor(M[istart, 4])
    n2 = np.ceil(M[islut, 4])
    a = M[istart, 4] - n1
    b = M[islut, 4] - np.floor(M[islut, 4])
    
    col1 = La1
    col2 = np.array((1 - a) * P[:, int(n1)] + a * P[:, int(n1)+1])
    if int(n1+1) == int(n2-1):
        col3 = np.array([P[:, int(n1+1)]])
    else:
        col3 = np.array(P[:, int(n1+1):int(n2)])
    col4 = np.array((1 - b) * P[:, int(n2 - 1)] + b * P[:, int(n2)])
    col5 = Lb2
    if len(col3) == 1:
        pts1 = np.concatenate([col1.reshape(-1,1), col2.reshape(-1,1), col3.reshape(-1,1), col4.reshape(-1,1), col5.reshape(-1,1)],axis=1)
    else:
        pts1 = np.concatenate([col1.reshape(-1,1), col2.reshape(-1,1), col3, col4.reshape(-1,1), col5.reshape(-1,1)],axis=1)

    if M[islut, 3] < M[istart, 3]:
        n3 = np.floor(M[islut, 3])
        n4 = np.ceil(M[istart, 3])
        a = M[islut, 3] - np.floor(M[islut, 3])
        b = M[istart, 3] - np.floor(M[istart, 3])
        
        col1 = ((1 - a) * P[:, int(n3)] + a * P[:, int(n3+1)]).reshape(-1,1)
        if int(n3+1) == int(n4-1):
            col2 = np.array([P[:, int(n3+1)]])
        else:
            col2 = (P[:, int(n3+1):int(n4)])
        col3 = ((1 - b) * P[:, int(n4 - 1)] + b * P[:, int(n4)]).reshape(-1,1)
        if len(col2) == 1:
            pts2 = np.concatenate((col1, col3), axis=1)
        else: 
            pts2 = np.concatenate((col1, col2, col3), axis=1)
    else:
        n3 = np.floor(M[istart, 3])
        n4 = np.ceil(M[islut, 3])
        a = M[istart, 3] - np.floor(M[istart, 3])
        b = M[islut, 3] - np.floor(M[islut, 3])

        col1 = ((1 - a) * P[:, int(n3)] + a * P[:, int(n3+1)]).reshape(-1,1)
        if int(n3+1) == int(n4-1):
            col2 = np.array([P[:, int(n3+1)]]).reshape(-1,1)
        else:
            col2 = (P[:, int(n3+1):int(n4)])
        col3 = ((1 - b) * P[:, int(n4 -  1)] + b * P[:, int(n4)]).reshape(-1,1)
        
        pts2 = np.concatenate((col1, col2, col3), axis=1)
        n = pts2.shape[1]
        pts2 = pts2[:, n::-1]
    ns = [n3, n4]
    n3 = np.min(ns)
    n4 = np.max(ns)
    pts = np.concatenate((pts1, pts2), axis=1)

    center = (np.sum(pts, axis=1) / pts.shape[1]).reshape(-1,1)
    pts = pts - np.tile(center, (1, pts.shape[1]))
    rdisk = np.max(np.sum(pts ** 2, axis=0)) ** 0.5
    pts = np.concatenate((pts, pts[:, 0].reshape(-1, 1)), axis=1)
    
    
    NbrTriangles = pts.shape[1] - 1
    P = P - np.tile(center, (1, P.shape[1]))
    Lindex = np.concatenate((np.arange(0, int(n1)), np.arange(int(n2), int(n3 - 1)), np.arange(int(n4), P.shape[1] - 1)))
    Lstart = P[:, Lindex]
    Lend = P[:, 1 + Lindex]

    Lmidt = np.sum(((Lstart + Lend) / 2) ** 2,axis=0) ** 0.5
    LineSegmentLength = np.sum((Lstart - Lend) ** 2,axis=0) ** 0.5
    distdiff = Lmidt - LineSegmentLength / 2
    ex = np.where(distdiff <= rdisk)[0]

    Lstart = Lstart[:, ex]
    Lend = Lend[:, ex]
    NbrL = Lstart.shape[1]
    NbrIntc = 0

    index = np.concatenate((np.arange(0, int(min(na, nb))), np.arange(int(max(na, nb))+1, int(min(na2, nb2)) - 1), np.arange(int(max(na2, nb2))+1, Pts1.shape[1] - 1)))
    Lstart1 = Pts1[:, index]
    Lend1 = Pts1[:, 1 + index]
    Lstart2 = Pts2[:, index]
    Lend2 = Pts2[:, 1 + index]

    LhomotopyCenter = 0.25 * (Lstart1 + Lend1 + Lstart2 + Lend2)
    LaCenter = 0.5 * (La1 + La2)
    LbCenter = 0.5 * (Lb1 + Lb2)
    La2Center = 0.5 * (La12 + La22)
    Lb2Center = 0.5 * (Lb12 + Lb22)
    radiusa = 0.5 * (np.sum((La1 - La2) ** 2)) ** 0.5
    radiusb = 0.5 * (np.sum((Lb1 - Lb2) ** 2)) ** 0.5
    radiusa2 = 0.5 * (np.sum((La12 - La22) ** 2)) ** 0.5
    radiusb2 = 0.5 * (np.sum((Lb12 - Lb22) ** 2)) ** 0.5

    rh1s = np.sum((Lstart1 - LhomotopyCenter) ** 2,axis=0)
    rh2s = np.sum((Lstart2 - LhomotopyCenter) ** 2,axis=0)
    rh1e = np.sum((Lend1 - LhomotopyCenter) ** 2,axis=0)
    rh2e = np.sum((Lend2 - LhomotopyCenter) ** 2,axis=0)
    radiushomotopy = np.amax(np.array([rh1s, rh2s, rh1e, rh2e]),axis=0) ** 0.5

    dista = np.sum((np.tile(LaCenter.reshape(-1,1), (1, len(index))) - LhomotopyCenter) ** 2,axis=0) ** 0.5
    distb = np.sum((np.tile(LbCenter.reshape(-1,1), (1, len(index))) - LhomotopyCenter) ** 2,axis=0) ** 0.5
    dista2 = np.sum((np.tile(La2Center.reshape(-1,1), (1, len(index))) - LhomotopyCenter) ** 2,axis=0) ** 0.5
    distb2 = np.sum((np.tile(Lb2Center.reshape(-1,1), (1, len(index))) - LhomotopyCenter) ** 2,axis=0) ** 0.5

    exa = np.where(dista <= radiusa + radiushomotopy)[0]
    exb = np.where(distb <= radiusb + radiushomotopy)[0]
    exa2 = np.where(dista2 <= radiusa2 + radiushomotopy)[0]
    exb2 = np.where(distb2 <= radiusb2 + radiushomotopy)[0]

    for ii in exa:
        tmp = itls.intersection_triangle_line_segment(Lstart1[:, ii], Lend1[:, ii], Lstart2[:, ii], La1, La2)[0] + itls.intersection_triangle_line_segment(Lend1[:, ii], Lstart2[:, ii], Lend2[:, ii], La1, La2)[0]
        NbrIntc += tmp
    for ii in exb:
        tmp = itls.intersection_triangle_line_segment(Lstart1[:, ii], Lend1[:, ii], Lstart2[:, ii], Lb1, Lb2)[0] + itls.intersection_triangle_line_segment(Lend1[:, ii], Lstart2[:, ii], Lend2[:, ii], Lb1, Lb2)[0]
        NbrIntc += tmp
    for ii in exa2:
        tmp = itls.intersection_triangle_line_segment(Lstart1[:, ii], Lend1[:, ii], Lstart2[:, ii], La12, La22)[0] + itls.intersection_triangle_line_segment(Lend1[:, ii], Lstart2[:, ii], Lend2[:, ii], La12, La22)[0]
        NbrIntc += tmp
    for ii in exb2:
        tmp = itls.intersection_triangle_line_segment(Lstart1[:, ii], Lend1[:, ii], Lstart2[:, ii], Lb12, Lb22)[0] + itls.intersection_triangle_line_segment(Lend1[:, ii], Lstart2[:, ii], Lend2[:, ii], Lb12, Lb22)[0]
        NbrIntc += tmp

    if NbrIntc > 0:
        ud = [0, 0]
        return ud

    dists = dp2l.d_points2line(pts, pts[:, 0], pts[:, pts1.shape[1]-1])
    ud = [np.sum(dists) * 2, lengRep]

    return ud

# M = np.array([
#     [54, 106, 52, 106.752476278650, 52.3523489677478, -1, 3.40061261331869, 0.570506487806912],
#     [16, 84, 68, 84.1318531700235, 68.5283480854625, -1, 0.279770334970397, 0.671405969257914],
#     [21, 89, 68, 89.0340661520450, 68.3484827537619, 1, 0.491691699666590, 0.627887069827490],
#     [4, 88, 84, 88.2008275408357, 84.6272145922403, -1, 0.0146748055529450, 0.654515355259122]
# ])

# M0 = np.array([
#     [54, 102, 48, 102.752476278650, 48.3523489677478, -1, 3.40061261331869, 0.570506487806912],
#     [16, 80, 64, 80.1318531700235, 64.5283480854625, -1, 0.279770334970397, 0.671405969257914],
#     [21, 85, 64, 85.0340661520450, 64.3484827537619, 1, 0.491691699666590, 0.627887069827490],
#     [4, 84, 80, 84.2008275408357, 80.6272145922403, -1, 0.0146748055529450, 0.654515355259122]
# ])

# M1 = np.array([
#     [54, 128, 74, 128.752476278650, 74.3523489677478, -1, 3.40061261331869, 0.570506487806912],
#     [16, 106, 90, 106.131853170023, 90.5283480854625, -1, 0.279770334970397, 0.671405969257914],
#     [21, 111, 90, 111.034066152045, 90.3484827537619, 1, 0.491691699666590, 0.627887069827490],
#     [4, 110, 106, 110.200827540836, 106.627214592240, -1, 0.0146748055529450, 0.654515355259122]
# ])

# i = 1
# j = 2
# makker = 3
# maxlen = 15
# P = np.loadtxt("/Users/agb/Desktop/Bachelor projekt/Python kode oversat/P.txt")
# P1 = np.loadtxt("/Users/agb/Desktop/Bachelor projekt/Python kode oversat/PP1.txt")

# tmp=IsContractableType2ReparametrizationParallel(M,M0,M1,i,j,P,P1,maxlen)

# i og makker (j) er 1 mindre end i Matlab grundet 0-indeksering!!!!!!
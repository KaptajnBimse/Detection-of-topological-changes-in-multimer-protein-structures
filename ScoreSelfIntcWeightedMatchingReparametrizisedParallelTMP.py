import numpy as np
from scipy import sparse
from scipy.interpolate import spmak
from scipy.interpolate import fnval
from scipy.optimize import linear_sum_assignment
from IsContractableType1ReparametrizationParallel import IsContractableType1ReparametrizationParallel
from IsContractableType2ReparametrizationParallel import IsContractableType2ReparametrizationParallel
from PriceEstEndContraction import PriceEstEndContraction
from distPP import distPP


def ScoreSelfIntcWeightedMatchingReparametrizisedParallelTMP(selfintc, selfintcu, selfintcv, selfintcs, len, P, P1, RePar1, RePar2, IsAligned, P1org, P2org, maxendcontraction, maxlen):
    # Ouptut order:
    #  1 number of local reidemeister type one moves performed
    #  2 total cost of local reidemeister type one moves performed
    #  3 total length of local reidemeister type one moves performed
    #  4 maximal cost of one local reidemeister type one moves performed
    #  5 maximal length of one local reidemeister type one moves performed

    #  6 number of local reidemeister type one moves performed
    #  7 total cost of local reidemeister type two moves performed
    #  8 total length of local reidemeister type two moves performed
    #  9 maximal cost of one local reidemeister type two moves performed
    # 10 maximal length of one local reidemeister type two moves performed
    # 11 morph length (L1) of enddeformations (pre and post morphs)
    # 12 morph length (L1) of the non-contracted part
    # 13 morph length (L1) of the ends during the original morph
    # 14 maximal depth of essencial self-intersections
    # 15 number of essencial self-intersections
    # 16 Sum of signs of the self-intersections
    # 17 Sum of signs of the self-intersections after local reidemeister moves

    # 18 Number of Essensial Aligned-Aligned self-intersections
    # 19 Number of Aligned-Aligned self-intersections
    # 20 Number of Essensial Aligned-Gap self-intersections
    # 21 Number of Aligned-Gap self-intersections
    # 22 Number of Essensial Gap-Gap self-intersections
    # 23 Number of Gap-Gap self-intersections

    # Re 18 to 23: A self-intersection between parameters (a,b) is defined as
    #         aligned-aligned if   IsAlignedSpline(a)+IsAlignedSpline(b) >= 1.5
    #         aligned-gap if 1.5 > IsAlignedSpline(a)+IsAlignedSpline(b) >  0.5
    #         gap-gap     if       IsAlignedSpline(a)+IsAlignedSpline(b) <  0.5

    C = sparse.find(sparse.tril(selfintc, 0))[2]
    d = sparse.find(sparse.tril(selfintcu, 0))[2]
    e = sparse.find(sparse.tril(selfintcv, 0))[2]
    f = sparse.find(sparse.tril(selfintcs, 0))[2]

    M = np.column_stack((C-d, C+e, C, C+d, C+e, C))
    M = np.column_stack((M, 0.001*(3.8**2*np.min([(M[:,4]-1)**2, (M[:,3]-M[:,4])**2/(4*np.pi), (len-M[:,3])**2], axis=0)), f))
    ud_M = M

    n1 = RePar1.shape[1]
    sp1 = spmak(np.arange(n1+2), RePar1[0])

    n2 = RePar2.shape[1]
    sp2 = spmak(np.arange(n2+2), RePar2[0])

    M0 = M.copy()
    M1 = M.copy()
    M0[:,1:4] = fnval(sp1, M[:,1:4])
    M1[:,1:4] = fnval(sp2, M[:,1:4])

    n3 = IsAligned.shape[1]
    IsAlignedSpline = spmak(np.arange(n3+2), IsAligned[0])

    Nbr = M.shape[0]
    NbrOmega1 = 0
    NbrOmega2 = 0
    cost1 = np.zeros(2)
    maxCost1 = np.zeros(2)
    cost2 = np.zeros(2)
    maxCost2 = np.zeros(2)
    sumsignraw = np.sum(M[:,5]) # sum of crossing changes
    NbrSelfIntc = M.shape[0]
    O1 = np.zeros((NbrSelfIntc, 2))
    for j in range(NbrSelfIntc):
        tmp = IsContractableType1ReparametrizationParallel(M, M0, M1, j, P, P1, maxlen)
        if tmp[0]:
            tmp[0] = np.min([tmp[0], PriceEstEndContraction(M[j,4]-1), PriceEstEndContraction(len-M[j,3])])
        else:
            enddist = np.min([M[j,4]-1, len-M[j,3]])
            if enddist < maxendcontraction:
                tmp = [PriceEstEndContraction(enddist), enddist*2]
        O1[j,:] = tmp

    paircount = 0
    O2 = np.zeros((Nbr*(Nbr-1)//2, 4))
    for i in range(Nbr-1):
        for j in range(i+1, Nbr):
            if M[i,5]+M[j,5] == 0: # have opposite signs
                if not (M[j,3] < M[i,4] or M[j,4] > M[i,3]):
                    tmp = IsContractableType2ReparametrizationParallel(M, M0, M1, i, j, P, P1, maxlen)
                    if tmp[0]:
                        paircount += 1
                        O2[paircount-1,:] = [i, j] + tmp

    O2 = O2[:paircount,:]

    epsilon = 0.5*(np.sum(O1[:,0]) + np.sum(O2[:,2]))**(-1)
    WVertex = epsilon*O1[:,0] + (O1[:,0] == 0)

    Wedge = -epsilon*O2[:,2] + WVertex[O2[:,0]] + WVertex[O2[:,1]]

    edgeData = np.column_stack((O2[:,0:2], Wedge))
    result = linear_sum_assignment(-edgeData)
    NbrEssensial = 0
    Essensials = []
    Nbr2 = len(result[0])

    for i in range(Nbr2):
        if result[0][i] > 0:
            if result[0][i] < i:
                edge = np.where((O2[:,0] == result[0][i]) & (O2[:,1] == i))[0][0]
                cost2var = O2[edge,2:4]
                cost2 += cost2var
                NbrOmega2 += 1
                maxCost2 = np.maximum(maxCost2, cost2var)
        else:
            if WVertex[i] < 1:
                cost1var = O1[i,:]
                cost1 += cost1var
                NbrOmega1 += 1
                maxCost1 = np.maximum(maxCost1, cost1var)
            else:
                NbrEssensial += 1
                Essensials.append(i)

    for i in range(Nbr2, Nbr):
        if WVertex[i] < 1:
            cost1var = O1[i,:]
            cost1 += cost1var
            NbrOmega1 += 1
            maxCost1 = np.maximum(maxCost1, cost1var)
        else:
            NbrEssensial += 1
            Essensials.append(i)

    slet = fnval(IsAlignedSpline, M[:,3]) + fnval(IsAlignedSpline, M[:,4])
    AlignedAligned = slet >= 1.5
    AlignedGap = (slet < 1.5) & (slet > 0.5)
    GapGap = slet <= 0.5
    NbrEssensialAlignedAligned = np.sum(AlignedAligned[Essensials])
    NbrAlignedAlignedTotal = np.sum(AlignedAligned)
    NbrEssensialAlignedGap = np.sum(AlignedGap[Essensials])
    NbrAlignedGapTotal = np.sum(AlignedGap)
    NbrEssensialGapGap = np.sum(GapGap[Essensials])
    NbrGapGapTotal = np.sum(GapGap)


    
    PerformEndDeformations = 0

    if PerformEndDeformations == 0:
        RMSsum = np.sum(np.sqrt(np.sum((P - P1) ** 2, axis=1)))
        if NbrEssensial == 0:
            ud = [NbrOmega1, cost1, maxCost1, NbrOmega2, cost2, maxCost2, 0, RMSsum, 0, 0, 0,
                  sumsignraw, 0, NbrEssensialAlignedAligned, NbrAlignedAlignedTotal,
                  NbrEssensialAlignedGap, NbrAlignedGapTotal, NbrEssensialGapGap,
                  NbrGapGapTotal]
            ud_essensials = np.zeros((0, 2))
        else:
            M = M[Essensials, :]
            ud_essensials = M[:, [1, 2]]
            antal = M.shape[0]
            ud = [NbrOmega1, cost1, maxCost1, NbrOmega2, cost2, maxCost2, 0,
                  RMSsum, 0, 0, antal, sumsignraw, np.sum(M[:, 5]),
                  NbrEssensialAlignedAligned, NbrAlignedAlignedTotal,
                  NbrEssensialAlignedGap, NbrAlignedGapTotal, NbrEssensialGapGap,
                  NbrGapGapTotal]
    else:
        if NbrEssensial == 0:
            RMSsum = np.sum(np.sqrt(np.sum((P - P1) ** 2, axis=1)))
            ud = [NbrOmega1, cost1, maxCost1, NbrOmega2, cost2, maxCost2, 0, RMSsum, 0, 0, 0,
                  sumsignraw, 0, NbrEssensialAlignedAligned, NbrAlignedAlignedTotal,
                  NbrEssensialAlignedGap, NbrAlignedGapTotal, NbrEssensialGapGap,
                  NbrGapGapTotal]
            ud_essensials = np.zeros((0, 2))
        # else:
        #     M = M[Essensials, :]
        #     a, b = ApproxCheapestEndContraction(M, len)
        #     ud_essensials = M[:, [1, 2]]
        #     enddeformations = 0
        #     P1orgWindow = P1org[RePar1[0, 0]:RePar1[0, -1], :]  # cutting the aligned windows
        #     P2orgWindow = P2org[RePar2[0, 0]:RePar2[0, -1], :]

        #     if a > 1:
        #         a1 = fnval(sp1, a)
        #         a2 = fnval(sp2, a)
        #         enddeformations += ContractStart(a1 - RePar1[0, 0] + 1, P1orgWindow.T) + \
        #                            ContractStart(a2 - RePar2[0, 0] + 1, P2orgWindow.T)
        #     if b < len:
        #         b1 = fnval(sp1, b)
        #         b2 = fnval(sp2, b)
        #         enddeformations += ContractEnd(b1 - RePar1[0, 0] + 1, P1orgWindow.T) + \
        #                            ContractEnd(b2 - RePar2[0, 0] + 1, P2orgWindow.T)

            # sa = a - np.floor(a)
            # sb = b - np.floor(b)
            # Pa = (1 - sa) * P[np.floor(a), :] + sa * P[np.ceil(a), :]
            # P1a = (1 - sa) * P1[np.floor(a), :] + sa * P1[np.ceil(a), :]
            # Pb = (1 - sb) * P[np.floor(b), :] + sb * P[np.ceil(b), :]
            # P1b = (1 - sb) * P1[np.floor(b), :] + sb * P1[np.ceil(b), :]
            # RMSsum = np.sum(np.sqrt(np.sum((P[np.ceil(a):np.floor(b), :] - P1[np.ceil(a):np.floor(b), :]) ** 2, axis=1)))
            # EndeFlytning = distPP(Pa.T, P1a.T) * np.floor(a) + distPP(Pb.T, P1b.T) * (len - np.floor(b))

            # antal = M.shape[0]
            # dybde = np.zeros((antal, 1))
            # for i in range(antal):
            #     dybde[i, 0] = np.sum((M[:, 5] <= M[i, 5]) * (M[i, 5] <= M[:, 4]))
            # ud = [NbrOmega1, cost1, maxCost1, NbrOmega2, cost2, maxCost2, enddeformations,
            #       RMSsum, EndeFlytning, np.max(dybde), antal, sumsignraw, np.sum(M[:, 6]),
            #       NbrEssensialAlignedAligned, NbrAlignedAlignedTotal, NbrEssensialAlignedGap,
            #       NbrAlignedGapTotal, NbrEssensialGapGap, NbrGapGapTotal]
            
    
    return ud, ud_essensials, ud_M




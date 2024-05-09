import numpy as np
from NEAMReparametrizationParallel import NEAMReparametrizationParallel
from ScoreSelfIntcWeightedMatchingReparametrizisedParallelTMP import ScoreSelfIntcWeightedMatchingReparametrizisedParallelTMP
from AlignmentMetaData import AlignmentMetaData
from SelfintersectionTransversal import SelfintersectionTransversal
from MakeSelfIntcFigureV3 import MakeSelfIntcFigureV3
import copy

def OverlapandSelfintersectParallelV3(P1Less4, P2Less4, RePar1Less4, RePar2Less4, IsAligned, P1org, P2org, NresAverage, options, False_lines, P1, P2, RePar1, RePar2, IsAligned_org, Insert_points_P1, Insert_points_P):
    Smoothning = options['Smoothning']
    AllowEndContractions = options['AllowEndContractions']
    AllMaxLengths = options['MaxLength']
    makefigure = options['MakeFigures']

    AlignmentMetaDataOut = AlignmentMetaData(RePar1Less4, RePar2Less4, IsAligned)

    n = len(P1Less4)
    m = len(P2Less4)
    if abs(n - m) > 0:
        print('Unequal sized protein structures intented superimposed')
        return

    bands = np.arange(1, 6)
    sumselfintc = np.zeros(len(bands))
    sumoverlap = np.zeros(len(bands))

    dPsq = (P1Less4 - P2Less4) ** 2  # working zone

    Dsqr = np.sum(dPsq, axis=1)
    Ds = np.sqrt(Dsqr)
    GDT_TS = (np.sum(Ds <= 1) + np.sum(Ds <= 2) + np.sum(Ds <= 4) + np.sum(Ds <= 8)) / (4 * n)
    d0sqr = (1.24 * (NresAverage - 15) ** (1.0 / 3.0) - 1.8) ** 2
    TM = np.sum(1.0 / (1.0 + Dsqr[IsAligned == 1] / d0sqr)) / NresAverage
    NbrAlignedXXX = np.sum(IsAligned == 1)
    rms1 = np.sum(Ds)
    rms2 = np.sqrt(np.sum(Dsqr) / n)

    rms1Aligned = np.sum(Ds[IsAligned == 1])
    rms2Aligned = np.sqrt(np.sum(Dsqr[IsAligned == 1]) / np.sum(IsAligned))

    overlap, _, _, _ = NEAMReparametrizationParallel(P1Less4, P2Less4, RePar1Less4, RePar2Less4, IsAligned, Smoothning)

    L1 = np.sqrt(np.sum((P1Less4[0:n - 1, :] - P1Less4[1:n, :]) ** 2, axis=1))
    L2 = np.sqrt(np.sum((P2Less4[0:n - 1, :] - P2Less4[1:n, :]) ** 2, axis=1))
    # histogram of L1 and L2
    #import matplotlib.pyplot as plt
    #plt.hist(L1, bins=300)
    #plt.hist(L2, bins=300)
    #plt.show()

    if Smoothning == 1:
        MaxL = 3.5
        MaxSum = 2.1
    else:
        MaxL = 4.0
        MaxSum = 2.5
    LmaxOK = np.maximum(L1, L2)
    if Smoothning == 1:  # compensating that the first two and the last two C alphas
        # not are changed by the smoothing operation 
        LmaxOK[0:2] = LmaxOK[0:2] - [0.5, 0.35]
        LmaxOK[-2:] = LmaxOK[-2:] - [0.35, 0.5]
    
    M = np.tile(LmaxOK, (n - 1, 1))
    M = np.maximum(M, M.T)

    selfintc = np.zeros((n-1, n-1))
    selfintcu = np.zeros((n-1, n-1))
    selfintcv = np.zeros((n-1, n-1))
    selfintcs = np.zeros((n-1, n-1))

    tmp = np.zeros((n-1, n-1, 4))
    tmp[:,:,0] = overlap[0:n-1, 0:n-1]
    tmp[:,:,1] = overlap[1:n, 0:n-1]
    tmp[:,:,2] = overlap[0:n-1, 1:n]
    tmp[:,:,3] = overlap[1:n, 1:n]
    Oav = np.sum(tmp, axis=2)

    a2, a1 = np.where(np.transpose(np.tril(Oav, -2) > MaxSum) + (np.tril(M, -2) > MaxL))
    #a2, a1 = np.where(np.transpose((np.tril(M, -2) > MaxL)))

    tjekliste = np.column_stack((a1, a2))
    
    tmp_num_check = tjekliste.shape[0]

    # Create a boolean mask where True indicates the numbers to keep
    mask = ~np.isin(tjekliste[:, 0], False_lines)

    # Apply the mask to filter out the rows
    tjekliste = tjekliste[mask]

    mask = ~np.isin(tjekliste[:, 1], False_lines)

    # Apply the mask to filter out the rows
    tjekliste = tjekliste[mask]

    PotSelfIntc = tjekliste.shape[0]
    print("Number to check: ", PotSelfIntc)
    print("Remove because false lines ", tmp_num_check - PotSelfIntc)

    Insert_points_P_tot = np.concatenate(list(Insert_points_P.values()), axis = 0)
    Insert_cumsum = np.cumsum(Insert_points_P_tot)
    IPP0_tjek = Insert_cumsum[tjekliste[:, 0]]
    IPP1_tjek = Insert_cumsum[tjekliste[:, 1]]
    
    P1_tot = np.concatenate(list(P1.values()), axis = 0)
    P2_tot = np.concatenate(list(P2.values()), axis = 0)

    selfintcI =[]
    selfintcJ =[]
    for k in range(tjekliste.shape[0]):
        i = int(tjekliste[k, 0] - IPP0_tjek[k])
        j = int(tjekliste[k, 1] - IPP1_tjek[k])
        # print(i,j)
        if (j == i+1) or (j == i-1):
            print("Removed because of adjacent chains")
        else:
            UdSelf = SelfintersectionTransversal(P1_tot[i:(i+2), :].T, P2_tot[i:(i+2), :].T, P1_tot[j:(j+2), :].T, P2_tot[j:(j+2), :].T)
            UdSelf = np.atleast_2d(UdSelf)
            selfintc[i, j] = UdSelf[0, 0]
            print(f"{k/(PotSelfIntc-1)*100:.2f}%")
            if UdSelf[0, 0] ** 2 == 1:
                selfintcu[i, j] = UdSelf[0, 1]
                selfintcv[i, j] = UdSelf[0, 2]
                selfintcs[i, j] = UdSelf[0, 3]
                selfintcI = np.append(selfintcI, i)
                selfintcJ = np.append(selfintcJ, j)
    print(len(np.where(selfintc)[0]))
    
    # plot 10 random lines that should intersect -----------------------------------

    # #generate 10 random numbers in the range 0-len(np.where(selfintc)[0])
    # random_numbers = np.random.randint(0, len(np.where(selfintc)[0]), 10)
    # import plotly.graph_objects as go

    # for i in random_numbers:
    #     line1 = np.where(selfintc)[0][i]
    #     line2 = np.where(selfintc)[1][i]
    #     fig = go.Figure()
    #     fig.add_trace(go.Scatter3d(x=((1-selfintcs[line1, line2])*P1_tot[line1:line1+2]+selfintcs[line1, line2]*P2_tot[line1:line1+2])[:,0].tolist(), 
    #                                y=((1-selfintcs[line1, line2])*P1_tot[line1:line1+2]+selfintcs[line1, line2]*P2_tot[line1:line1+2])[:,1].tolist(), 
    #                                z=((1-selfintcs[line1, line2])*P1_tot[line1:line1+2]+selfintcs[line1, line2]*P2_tot[line1:line1+2])[:,2].tolist(),
    #                                mode='lines', line=dict(width=9,color = 'yellow'), name='line 1'))
        
    #     fig.add_trace(go.Scatter3d(x=((1-(selfintcs[line1, line2]+0.1))*P1_tot[line1:line1+2]+(selfintcs[line1, line2]+0.1)*P2_tot[line1:line1+2])[:,0].tolist(), 
    #                                y=((1-(selfintcs[line1, line2]+0.1))*P1_tot[line1:line1+2]+(selfintcs[line1, line2]+0.1)*P2_tot[line1:line1+2])[:,1].tolist(), 
    #                                z=((1-(selfintcs[line1, line2]+0.1))*P1_tot[line1:line1+2]+(selfintcs[line1, line2]+0.1)*P2_tot[line1:line1+2])[:,2].tolist(), 
    #                                mode='lines', line=dict(width=9,color = 'blue'), name='line 1'))

    #     fig.add_trace(go.Scatter3d(x=((1-(selfintcs[line1, line2]-0.1))*P1_tot[line1:line1+2]+(selfintcs[line1, line2]-0.1)*P2_tot[line1:line1+2])[:,0].tolist(), 
    #                                y=((1-(selfintcs[line1, line2]-0.1))*P1_tot[line1:line1+2]+(selfintcs[line1, line2]-0.1)*P2_tot[line1:line1+2])[:,1].tolist(), 
    #                                z=((1-(selfintcs[line1, line2]-0.1))*P1_tot[line1:line1+2]+(selfintcs[line1, line2]-0.1)*P2_tot[line1:line1+2])[:,2].tolist(), 
    #                                mode='lines', line=dict(width=9,color = 'blue'), name='line 1'))

    #     fig.add_trace(go.Scatter3d(x=((1-selfintcs[line1, line2])*P1_tot[line2:line2+2]+selfintcs[line1, line2]*P2_tot[line2:line2+2])[:,0].tolist(), 
    #                                y=((1-selfintcs[line1, line2])*P1_tot[line2:line2+2]+selfintcs[line1, line2]*P2_tot[line2:line2+2])[:,1].tolist(), 
    #                                z=((1-selfintcs[line1, line2])*P1_tot[line2:line2+2]+selfintcs[line1, line2]*P2_tot[line2:line2+2])[:,2].tolist(), 
    #                                mode='lines', line=dict(width=9,color = 'yellow'), name='line 1'))
        
    #     fig.add_trace(go.Scatter3d(x=((1-(selfintcs[line1, line2]+0.1))*P1_tot[line2:line2+2]+(selfintcs[line1, line2]+0.1)*P2_tot[line2:line2+2])[:,0].tolist(), 
    #                                y=((1-(selfintcs[line1, line2]+0.1))*P1_tot[line2:line2+2]+(selfintcs[line1, line2]+0.1)*P2_tot[line2:line2+2])[:,1].tolist(), 
    #                                z=((1-(selfintcs[line1, line2]+0.1))*P1_tot[line2:line2+2]+(selfintcs[line1, line2]+0.1)*P2_tot[line2:line2+2])[:,2].tolist(), 
    #                                mode='lines', line=dict(width=9,color = 'red'), name='line 1'))

    #     fig.add_trace(go.Scatter3d(x=((1-(selfintcs[line1, line2]-0.1))*P1_tot[line2:line2+2]+(selfintcs[line1, line2]-0.1)*P2_tot[line2:line2+2])[:,0].tolist(), 
    #                                y=((1-(selfintcs[line1, line2]-0.1))*P1_tot[line2:line2+2]+(selfintcs[line1, line2]-0.1)*P2_tot[line2:line2+2])[:,1].tolist(), 
    #                                z=((1-(selfintcs[line1, line2]-0.1))*P1_tot[line2:line2+2]+(selfintcs[line1, line2]-0.1)*P2_tot[line2:line2+2])[:,2].tolist(), 
    #                                mode='lines', line=dict(width=9,color = 'red'), name='line 1'))

    #     fig.show()

    # --------------------------------------------------------------------------------

    for j in range(len(bands)):
        sumoverlap[j] = np.sum(np.tril(overlap, -bands[j]))
        sumselfintc[j] = np.sum(np.tril(np.abs(selfintc), -bands[j]))

    Maxs = AllMaxLengths
    Outs = []
    Udessentials = np.zeros((1,2))
    Intersecting_chain_number_i = np.zeros((1))
    Intersecting_chain_number_j = np.zeros((1))
    
    # Find the i,j in tjekliste where selfintersection is present (not intersection between chains)
    intersect_index_i = np.where(selfintc)[0]
    intersect_index_j = np.where(selfintc)[1]
    
    chain_change = np.zeros(len(P1)+1)
    start = -1
    for i,chain in zip(range(len(P1)), P1.keys()):
        chain_change[i+1] = len(P1[chain])+start
        start = chain_change[i+1]
    chain_change2 = copy.deepcopy(chain_change)
    chain_change = np.delete(chain_change, -1) # Remove the last element because for loop should not include it
    
    selfintersect = {}
    non_selfintersect = {}
    for c in reversed(chain_change):
        cc = int(np.where(chain_change == c)[0])
        chain = list(P1.keys())[cc]
        selfintersect[chain] = []
        non_selfintersect[chain] = []
        for i, k in enumerate(reversed(range(intersect_index_i.shape[0]))):
            if intersect_index_i[k] > c and intersect_index_j[k] > c:
                #selfintersect[chain][i] = ([intersect_index_i[k], intersect_index_j[k]]) # If you want to know indices
                selfintersect[chain].append(([intersect_index_i[k], intersect_index_j[k]]))
                intersect_index_i = np.delete(intersect_index_i, k)
                intersect_index_j = np.delete(intersect_index_j, k)
                continue
            if (intersect_index_i[k] < c) ^ (intersect_index_j[k] < c):
                #non_selfintersect[chain][i] = ([intersect_index_i[k], intersect_index_j[k]])
                non_selfintersect[chain].append(([intersect_index_i[k], intersect_index_j[k]]))
                intersect_index_i = np.delete(intersect_index_i, k)
                intersect_index_j = np.delete(intersect_index_j, k)
    
    selfintersect_tot = []
    non_selfintersect_tot = []
    
    
    
    
    for i in list(selfintersect.keys()):
        selfintersect_tot.extend(selfintersect[i])
        non_selfintersect_tot.extend(non_selfintersect[i])


    # Assuming selfintc and selfintersect[chain] are defined
    # Reset all elements in selfintc to 0
    for i in range(np.array([Maxs]).shape[0]):
        if AllowEndContractions == 1:
            maxendcontraction = Maxs[i] / 2
        else:
            maxendcontraction = 0
        
        for i in range(len(P1.keys())):
            for j in range(i,len(P1.keys())):
                new_selfintc = np.zeros((selfintc.shape[0], selfintc.shape[1]))
                new_selfintcu = np.zeros((selfintc.shape[0], selfintc.shape[1]))
                new_selfintcv = np.zeros((selfintc.shape[0], selfintc.shape[1]))
                new_selfintcs = np.zeros((selfintc.shape[0], selfintc.shape[1]))
                startj = int(chain_change2[i])
                starti = int(chain_change2[j])
                slutj = int(chain_change2[i+1])
                sluti = int(chain_change2[j+1])
                new_selfintc[starti:sluti, startj:slutj] = selfintc[starti:sluti, startj:slutj]
                new_selfintcu[starti:sluti, startj:slutj] = selfintcu[starti:sluti, startj:slutj]
                new_selfintcv[starti:sluti, startj:slutj] = selfintcv[starti:sluti, startj:slutj]
                new_selfintcs[starti:sluti, startj:slutj] = selfintcs[starti:sluti, startj:slutj]
                if np.where(new_selfintc)[0].shape[0] != 0:
                    print(i,j)
                    tmp, Essensials, Mselfintc = ScoreSelfIntcWeightedMatchingReparametrizisedParallelTMP(new_selfintc, new_selfintcu, new_selfintcv, new_selfintcs, n, P1_tot, P2_tot, RePar1, RePar2, IsAligned, i, j, maxendcontraction, Maxs)
                    if len(Essensials) != 0:
                        Udessentials = np.vstack((Udessentials, Essensials))
                        Intersecting_chain_number_i = np.hstack((Intersecting_chain_number_i, np.ones(Essensials.shape[0])*i))
                        Intersecting_chain_number_j = np.hstack((Intersecting_chain_number_j, np.ones(Essensials.shape[0])*j))

                    Outs.append(tmp)
    
    """"
    for i in range(Maxs):
        if AllowEndContractions == 1:
            maxendcontraction = Maxs[i] / 2
        else:
            maxendcontraction = 0
        
        tmp, Essensials, Mselfintc = ScoreSelfIntcWeightedMatchingReparametrizisedParallelTMP(selfintc, selfintcu, selfintcv, selfintcs, n, P1_tot, P2_tot, RePar1, RePar2, IsAligned, P1org, P2org, maxendcontraction, Maxs)
        Outs.append(tmp)
    """
    Udessentials = Udessentials[1:,:]
    if makefigure == 1:
        MakeSelfIntcFigureV3(P1_tot, P2_tot, selfintc, overlap, Udessentials, RePar1, RePar2, options, chain_change,Intersecting_chain_number_i ,Intersecting_chain_number_j)

    ud = [Outs, rms1, rms1Aligned, rms2, rms2Aligned, GDT_TS, TM, sumoverlap, PotSelfIntc, sumselfintc, AlignmentMetaDataOut]

# P1 = np.loadtxt('Monomer/Test txt/TestEssential/P1.txt')
# P2 = np.loadtxt('Monomer/Test txt/TestEssential/P2.txt')
# RePar1 = np.loadtxt('Monomer/Test txt/TestEssential/RePar1.txt')
# RePar2 = np.loadtxt('Monomer/Test txt/TestEssential/RePar2.txt')
# IsAligned = np.loadtxt('Monomer/Test txt/TestEssential/IsAligned.txt')
# P1org = np.loadtxt('Monomer/Test txt/TestEssential/P1org.txt')
# P2org = np.loadtxt('Monomer/Test txt/TestEssential/P2org.txt')
# NresAverage = np.loadtxt('Monomer/Test txt/TestEssential/NresAverage.txt')



# P1 = np.loadtxt('Monomer/Test txt/Omega2a_b/P1.txt')
# P2 = np.loadtxt('Monomer/Test txt/Omega2a_b/P2.txt')
# RePar1 = np.loadtxt('Monomer/Test txt/Omega2a_b/RePar1.txt')
# RePar2 = np.loadtxt('Monomer/Test txt/Omega2a_b/RePar2.txt')
# IsAligned = np.loadtxt('Monomer/Test txt/Omega2a_b/IsAligned.txt')
# P1org = np.loadtxt('Monomer/Test txt/Omega2a_b/P1org.txt')
# P2org = np.loadtxt('Monomer/Test txt/Omega2a_b/P2org.txt')
# NresAverage = np.loadtxt('Monomer/Test txt/Omega2a_b/NresAverage.txt')


# options_fig = {
#     'MaxLength': 15,
#     'dmax': 10,
#     'Smoothning': 0,
#     'AllowEndContractions': 0,
#     'MakeFigures': 1,
#     'MakeAlignmentSeedFigure': 0,
#     'MakeFiguresInLastItteration': 1,
#     'MakeLocalPlotsOfEssensials': 1,
#     'SelfIntcFigCutSize': 10,
#     'PrintOut': 0,
#     'additionalRMSD': 0,
#     'alignmentsmoothing': 0,
#     'alignmentsmoothingwidth': 3,
#     'AdaptiveSubset': 1,
#     'MaxNbrAlignmentSeeds': 7,
#     'MaxSeedOverlap': 0.5000,
#     'MinSeedLength': 40,
#     'OverlapWeight': 4,
#     'MaxIter': 20,
#     'MaxWindowMisalignment': 1,
#     'MaxMisAlignment': 0.0150,
#     'MinimalAlignmentLength': 30,
#     'FileName1': 'file1.pdb',
#     'FileName2': 'file2.pdb',
#     'StructureSequenceWeight': 1.5608,
#     'SeqenceMisAlignmentPenalty': [7.2200  ,  2.1660], 
#     'TrimSeqenceAlignment': 0,
#     'SequenceAlignmentExtension': 1,
#     'InitialAlignmentExactPairs': 1
#     }

# OverlapandSelfintersectParallelV3(P1, P2, RePar1, RePar2, IsAligned, P1org, P2org, NresAverage, options_fig)
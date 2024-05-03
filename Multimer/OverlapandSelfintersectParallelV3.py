import numpy as np
from NEAMReparametrizationParallel import NEAMReparametrizationParallel
from ScoreSelfIntcWeightedMatchingReparametrizisedParallelTMP import ScoreSelfIntcWeightedMatchingReparametrizisedParallelTMP
from AlignmentMetaData import AlignmentMetaData
from SelfintersectionTransversal import SelfintersectionTransversal
from MakeSelfIntcFigureV3 import MakeSelfIntcFigureV3

def OverlapandSelfintersectParallelV3(P1, P2, RePar1, RePar2, IsAligned, P1org, P2org, NresAverage, options, False_lines, P1Less4, P2Less4, RePar1Less4, ReParLess4):
    Smoothning = options['Smoothning']
    AllowEndContractions = options['AllowEndContractions']
    AllMaxLengths = options['MaxLength']
    makefigure = options['MakeFigures']

    AlignmentMetaDataOut = AlignmentMetaData(RePar1, RePar2, IsAligned)

    n = len(P1)
    m = len(P2)
    if abs(n - m) > 0:
        print('Unequal sized protein structures intented superimposed')
        return

    bands = np.arange(1, 6)
    sumselfintc = np.zeros(len(bands))
    sumoverlap = np.zeros(len(bands))

    dPsq = (P1 - P2) ** 2  # working zone

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

    overlap, _, _, _ = NEAMReparametrizationParallel(P1, P2, RePar1, RePar2, IsAligned, Smoothning)

    L1 = np.sqrt(np.sum((P1[0:n - 1, :] - P1[1:n, :]) ** 2, axis=1))
    L2 = np.sqrt(np.sum((P2[0:n - 1, :] - P2[1:n, :]) ** 2, axis=1))
    # histogram of L1 and L2
    import matplotlib.pyplot as plt
    plt.hist(L1, bins=300)
    plt.hist(L2, bins=300)
    plt.show()

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




    for k in range(tjekliste.shape[0]):
        i = tjekliste[k, 0]
        j = tjekliste[k, 1]
        UdSelf = SelfintersectionTransversal(P1[i:(i+2), :].T, P2[i:(i+2), :].T, P1[j:(j+2), :].T, P2[j:(j+2), :].T)
        UdSelf = np.atleast_2d(UdSelf)
        selfintc[i, j] = UdSelf[0, 0]
        print(f"{k/PotSelfIntc*100:.2f}%")
        if UdSelf[0, 0] ** 2 == 1:
            selfintcu[i, j] = UdSelf[0, 1]
            selfintcv[i, j] = UdSelf[0, 2]
            selfintcs[i, j] = UdSelf[0, 3]

    for j in range(len(bands)):
        sumoverlap[j] = np.sum(np.tril(overlap, -bands[j]))
        sumselfintc[j] = np.sum(np.tril(np.abs(selfintc), -bands[j]))

    Maxs = AllMaxLengths
    Outs = []
    
    for i in range(Maxs):
        if AllowEndContractions == 1:
            maxendcontraction = Maxs[i] / 2
        else:
            maxendcontraction = 0
        tmp, Essensials, Mselfintc = ScoreSelfIntcWeightedMatchingReparametrizisedParallelTMP(selfintc, selfintcu, selfintcv, selfintcs, n, P1, P2, RePar1, RePar2, IsAligned, P1org, P2org, maxendcontraction, Maxs)
        Outs.append(tmp)

    if makefigure == 1:
        MakeSelfIntcFigureV3(P1, P2, selfintc, overlap, Essensials, RePar1, RePar2, options)

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
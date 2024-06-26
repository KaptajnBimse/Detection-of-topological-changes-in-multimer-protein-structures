import numpy as np
from NEAMReparametrizationParallel import NEAMReparametrizationParallel
from ScoreSelfIntcWeightedMatchingReparametrizisedParallelTMP import ScoreSelfIntcWeightedMatchingReparametrizisedParallelTMP
from AlignmentMetaData import AlignmentMetaData
from SelfintersectionTransversal import SelfintersectionTransversal
from MakeSelfIntcFigureV3 import MakeSelfIntcFigureV3

def OverlapandSelfintersectParallelV3(P1_multi, P2_multi, RePar1_multi, RePar2_multi, IsAligned_multi, P1org, P2org, NresAverage_multi, options):
    Smoothning = options['Smoothning']
    AllowEndContractions = options['AllowEndContractions']
    AllMaxLengths = options['MaxLength']
    makefigure = options['MakeFigures']


    # AlignmentMetaDataOut = AlignmentMetaData(RePar1, RePar2, IsAligned)
    rms1 = {}
    rms2 = {}
    GDT_TS = {}
    TM = {}
    rms1Aligned = {}
    rms2Aligned = {}
    #overlap = {}
    #Essentials = {}
    #Outs = {}

    for i in range(len(P1_multi)):
        chain = list(P1_multi.keys())[i]
        P1 = P1_multi[chain]
        P2 = P2_multi[chain]
        RePar1 = RePar1_multi[chain]
        RePar2 = RePar2_multi[chain]
        IsAligned = IsAligned_multi[chain]
        NresAverage = NresAverage_multi[chain]

        n = len(P1)
        m = len(P2)
        if abs(n - m) > 0:
            print('Unequal sized protein structures intented superimposed')
            return



        dPsq = (P1 - P2) ** 2  # working zone

        Dsqr = np.sum(dPsq, axis=1)
        Ds = np.sqrt(Dsqr)
        GDT_TS[chain] = (np.sum(Ds <= 1) + np.sum(Ds <= 2) + np.sum(Ds <= 4) + np.sum(Ds <= 8)) / (4 * n)
        d0sqr = (1.24 * (NresAverage - 15) ** (1.0 / 3.0) - 1.8) ** 2
        TM[chain] = np.sum(1.0 / (1.0 + Dsqr[IsAligned == 1] / d0sqr)) / NresAverage
        NbrAlignedXXX = np.sum(IsAligned == 1)
        rms1[chain] = np.sum(Ds)
        rms2[chain] = np.sqrt(np.sum(Dsqr) / n)

        rms1Aligned[chain] = np.sum(Ds[IsAligned == 1])
        rms2Aligned[chain] = np.sqrt(np.sum(Dsqr[IsAligned == 1]) / np.sum(IsAligned))

    P1_tot = np.concatenate(list(P1_multi.values()), axis = 0)
    P2_tot = np.concatenate(list(P2_multi.values()), axis = 0)
    
    index1 = 0
    index2 = 0
    RePar1_tot = []
    RePar2_tot = []

    for i in list(RePar1_multi.keys()):
        RePar1_tot.extend(RePar1_multi[i]+np.ones(len(RePar1_multi[i]))*index1)
        index1 += len(RePar1_multi[i])
        RePar2_tot.extend(RePar2_multi[i]+np.ones(len(RePar2_multi[i]))*index2)
        index2 += len(RePar2_multi[i])

    IsAligned_tot = np.ones(len(RePar1_tot))
        
    # Save the array to a text file
    # np.savetxt('P1_tot.txt', P1_tot)
    # np.savetxt('P2_tot.txt', P2_tot)
    # np.savetxt('RePar1_tot.txt', RePar1_tot)
    # np.savetxt('RePar2_tot.txt', RePar2_tot)
    # np.savetxt('IsAligned_tot.txt', IsAligned_tot)

    overlap_tot, _, _, _ = NEAMReparametrizationParallel(P1_tot, P2_tot, RePar1_tot, RePar2_tot, IsAligned_tot, Smoothning)
    n = len(P1_tot)
    m = len(P2_tot)

    import plotly.express as px

    fig = px.imshow(overlap_tot, labels=dict(x="Residue 1", y="Residue 2", color="Overlap"), color_continuous_scale='Greys')
    fig.show()

    L1 = np.sqrt(np.sum((P1_tot[0:n - 1, :] - P1_tot[1:n, :]) ** 2, axis=1))
    L2 = np.sqrt(np.sum((P2_tot[0:n - 1, :] - P2_tot[1:n, :]) ** 2, axis=1))

    if Smoothning == 1:
        MaxL = 3.5
        MaxSum = 2.1
    else:
        MaxL = 4.0
        MaxSum = 2.5
    # ------------------------------------------------

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
    tmp[:,:,0] = overlap_tot[0:n-1, 0:n-1]
    tmp[:,:,1] = overlap_tot[1:n, 0:n-1]
    tmp[:,:,2] = overlap_tot[0:n-1, 1:n]
    tmp[:,:,3] = overlap_tot[1:n, 1:n]
    Oav = np.sum(tmp, axis=2)

    a2, a1 = np.where(np.transpose(np.tril(Oav, -2) > MaxSum) + (np.tril(M, -2) > MaxL))
    tjekliste = np.column_stack((a1, a2))
    

    start = len(P1_multi[list(P1_multi.keys())[0]])
    remove_nums = np.zeros(len(P1_multi)-1)
    for i in range(len(P1_multi)-1):
        remove_nums[i] = start-1
        start += len(P1_multi[list(P1_multi.keys())[i+1]])

    # remove_nums[-1] = 2856
    # Create a boolean mask where True indicates the numbers to keep
    mask = ~np.isin(tjekliste[:, 0], remove_nums)

    # Apply the mask to filter out the rows
    tjekliste = tjekliste[mask]

    mask = ~np.isin(tjekliste[:, 1], remove_nums)

    # Apply the mask to filter out the rows
    tjekliste = tjekliste[mask]

    PotSelfIntc = tjekliste.shape[0]

    for k in range(tjekliste.shape[0]):
        i = tjekliste[k, 0]
        j = tjekliste[k, 1]
        UdSelf = SelfintersectionTransversal(P1_tot[i:(i+2), :].T, P2_tot[i:(i+2), :].T, P1_tot[j:(j+2), :].T, P2_tot[j:(j+2), :].T)
        UdSelf = np.atleast_2d(UdSelf)
        selfintc[i, j] = UdSelf[0, 0]
        print(f"{k/PotSelfIntc*100:.2f}%")
        if UdSelf[0, 0] ** 2 == 1:
            selfintcu[i, j] = UdSelf[0, 1]
            selfintcv[i, j] = UdSelf[0, 2]
            selfintcs[i, j] = UdSelf[0, 3]

    bands = np.arange(1, 6)
    sumselfintc = np.zeros(len(bands))
    sumoverlap = np.zeros(len(bands))

    for j in range(len(bands)):
        sumoverlap[j] = np.sum(np.tril(overlap_tot, -bands[j]))
        sumselfintc[j] = np.sum(np.tril(np.abs(selfintc), -bands[j]))

    Maxs = AllMaxLengths
    Outs = []
    for i in range(Maxs):
        if AllowEndContractions == 1:
            maxendcontraction = Maxs[i] / 2
        else:
            maxendcontraction = 0
        tmp, Essensials, Mselfintc = ScoreSelfIntcWeightedMatchingReparametrizisedParallelTMP(selfintc, selfintcu, selfintcv, selfintcs, n, P1_tot, P2_tot, RePar1_tot, RePar2_tot, IsAligned_tot, P1org, P2org, maxendcontraction, Maxs)
        Outs.append(tmp)

    # if makefigure == 1:
    #     MakeSelfIntcFigureV3(P1, P2, selfintc, overlap, Essensials, RePar1, RePar2, options)

    ud = [Outs, rms1, rms1Aligned, rms2, rms2Aligned, GDT_TS, TM, sumoverlap, PotSelfIntc, sumselfintc]

    # ------------------------------------------------
    """for i in range(len(P1_multi)):
        chain = list(P1_multi.keys())[i]
        P1 = P1_multi[chain]
        P2 = P2_multi[chain]
        RePar1 = RePar1_multi[chain]
        RePar2 = RePar2_multi[chain]
        n = len(P1)
             

        L1 = np.sqrt(np.sum((P1[0:n - 1, :] - P1[1:n, :]) ** 2, axis=1))
        L2 = np.sqrt(np.sum((P2[0:n - 1, :] - P2[1:n, :]) ** 2, axis=1))

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

        start_inner = 0
        for j in range(i,len(P1_multi)):
            chain2 = list(P1_multi.keys())[j]
            P1_2 = P1_multi[chain2]
            P2_2 = P2_multi[chain2]
            m = len(P2_2)
            bands[chain+chain2] = np.arange(1, 6)
            sumselfintc[chain+chain2] = np.zeros(len(bands[chain+chain2]))
            sumoverlap[chain+chain2] = np.zeros(len(bands[chain+chain2]))

            overlap_tmp = overlap_tot[start_outer:start_outer+n, start_inner:start_inner+m]

            tmp = np.zeros((n-1, n-1, 4))
            tmp[:,:,0] = overlap_tmp[0:n-1, 0:n-1]
            tmp[:,:,1] = overlap_tmp[1:n, 0:n-1]
            tmp[:,:,2] = overlap_tmp[0:n-1, 1:n]
            tmp[:,:,3] = overlap_tmp[1:n, 1:n]
            Oav = np.sum(tmp, axis=2)

            a2, a1 = np.where(np.transpose(np.tril(Oav, -2) > MaxSum) + (np.tril(M, -2) > MaxL))
            tjekliste = np.column_stack((a1, a2))
            PotSelfIntc = tjekliste.shape[0]
            start_inner += len(P1_2)

            for k in range(tjekliste.shape[0]):
                i = tjekliste[k, 0]
                j = tjekliste[k, 1]
                UdSelf = SelfintersectionTransversal(P1[i:(i+2), :].T, P2[i:(i+2), :].T, P1[j:(j+2), :].T, P2[j:(j+2), :].T)
                UdSelf = np.atleast_2d(UdSelf)
                selfintc[i, j] = UdSelf[0, 0]
                if UdSelf[0, 0] ** 2 == 1:
                    selfintcu[i, j] = UdSelf[0, 1]
                    selfintcv[i, j] = UdSelf[0, 2]
                    selfintcs[i, j] = UdSelf[0, 3]

            for j in range(len(bands)):
                print("hej")
                sumoverlap[chain+chain2][j] = np.sum(np.tril(overlap_tmp, -bands[chain+chain2][j]))
                sumselfintc[chain+chain2][j] = np.sum(np.tril(np.abs(selfintc), -bands[chain+chain2][j]))

            Maxs = AllMaxLengths
            Outs[chain+chain2] = []
            for i in range(Maxs):
                if AllowEndContractions == 1:
                    maxendcontraction = Maxs[i] / 2
                else:
                    maxendcontraction = 0
                tmp, Essentials[chain+chain2], Mselfintc = ScoreSelfIntcWeightedMatchingReparametrizisedParallelTMP(selfintc, selfintcu, selfintcv, selfintcs, n, P1, P2, RePar1, RePar2, IsAligned, P1org, P2org, maxendcontraction, Maxs)
                Outs[chain+chain2].append(tmp)

    if makefigure == 1:
        MakeSelfIntcFigureV3(P1, P2, selfintc, overlap, Essentials, RePar1, RePar2, options)

    ud = [Outs, rms1, rms1Aligned, rms2, rms2Aligned, GDT_TS, TM, sumoverlap, PotSelfIntc, sumselfintc]"""

#P1 = np.loadtxt('Test txt/TestEssential/P1.txt')
#P2 = np.loadtxt('Test txt/TestEssential/P2.txt')
#RePar1 = np.loadtxt('Test txt/TestEssential/RePar1.txt')
#RePar2 = np.loadtxt('Test txt/TestEssential/RePar2.txt')
#IsAligned = np.loadtxt('Test txt/TestEssential/IsAligned.txt')
#P1org = np.loadtxt('Test txt/TestEssential/P1org.txt')
#P2org = np.loadtxt('Test txt/TestEssential/P2org.txt')
#NresAverage = np.loadtxt('Test txt/TestEssential/NresAverage.txt')



# P1 = np.loadtxt('Monomer/Test txt/Omega2a_b/P1.txt')
# P2 = np.loadtxt('Monomer/Test txt/Omega2a_b/P2.txt')
# RePar1 = np.loadtxt('Monomer/Test txt/Omega2a_b/RePar1.txt')
# RePar2 = np.loadtxt('Monomer/Test txt/Omega2a_b/RePar2.txt')
# IsAligned = np.loadtxt('Monomer/Test txt/Omega2a_b/IsAligned.txt')
# P1org = np.loadtxt('Monomer/Test txt/Omega2a_b/P1org.txt')
# P2org = np.loadtxt('Monomer/Test txt/Omega2a_b/P2org.txt')
# NresAverage = np.loadtxt('Monomer/Test txt/Omega2a_b/NresAverage.txt')


"""#options_fig = {
    'MaxLength': 15,
    'dmax': 10,
    'Smoothning': 0,
    'AllowEndContractions': 0,
    'MakeFigures': 1,
    'MakeAlignmentSeedFigure': 0,
    'MakeFiguresInLastItteration': 1,
    'MakeLocalPlotsOfEssensials': 1,
    'SelfIntcFigCutSize': 10,
    'PrintOut': 0,
    'additionalRMSD': 0,
    'alignmentsmoothing': 0,
    'alignmentsmoothingwidth': 3,
    'AdaptiveSubset': 1,
    'MaxNbrAlignmentSeeds': 7,
    'MaxSeedOverlap': 0.5000,
    'MinSeedLength': 40,
    'OverlapWeight': 4,
    'MaxIter': 20,
    'MaxWindowMisalignment': 1,
    'MaxMisAlignment': 0.0150,
    'MinimalAlignmentLength': 30,
    'FileName1': 'file1.pdb',
    'FileName2': 'file2.pdb',
    'StructureSequenceWeight': 1.5608,
    'SeqenceMisAlignmentPenalty': [7.2200  ,  2.1660], 
    'TrimSeqenceAlignment': 0,
    'SequenceAlignmentExtension': 1,
    'InitialAlignmentExactPairs': 1
    }"""

#OverlapandSelfintersectParallelV3(P1, P2, RePar1, RePar2, IsAligned, P1org, P2org, NresAverage, options_fig)
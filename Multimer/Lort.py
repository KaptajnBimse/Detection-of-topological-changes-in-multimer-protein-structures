import numpy as np
from Structural_AlignmentV2 import structural_alignment
from OverlapandSelfintersectParallelV3 import OverlapandSelfintersectParallelV3

#pdb_file1 = "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB//CRUA_hexamer_positive.pdb"
#pdb_file2 = "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRU1_hexamer_negative.pdb"

pdb_file1 = "/Users/agb/Desktop/Bachelor projekt/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRUA_hexamer_positive.pdb"
pdb_file2 = "/Users/agb/Desktop/Bachelor projekt/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRU1_hexamer_negative.pdb"

P1, P2, RePar1, RePar2, IsAligned, NresAverage, P1Less4, P2Less4, RePar1Less4, RePar2Less4, Insert_points_P1, Insert_points_P =  structural_alignment(pdb_file1, pdb_file2, makefigure = 0)
options = {'Smoothning': 0, 'AllowEndContractions': 0, 'MaxLength': 15, 'MakeFigures': 1}
P1org = 0
P2org = 0

P1_tot = np.concatenate(list(P1.values()), axis = 0)
P2_tot = np.concatenate(list(P2.values()), axis = 0)
P1Less4_tot = np.concatenate(list(P1Less4.values()), axis = 0)
P2Less4_tot = np.concatenate(list(P2Less4.values()), axis = 0)


index1 = 0
index2 = 0
index3 = 0
index4 = 0
RePar1_tot = []
RePar2_tot = []
RePar1Less4_tot = []
RePar2Less4_tot = []

for i in list(RePar2.keys()):
    RePar1_tot.extend(RePar1[i]+np.ones(len(RePar1[i]))*index1)
    index1 += RePar1[i][-1]+1
    RePar2_tot.extend(RePar2[i]+np.ones(len(RePar2[i]))*index2)
    index2 += RePar2[i][-1]+1
    RePar1Less4_tot.extend(RePar1Less4[i]+np.ones(len(RePar1Less4[i]))*index3)
    index3 += RePar1Less4[i][-1]+1
    RePar2Less4_tot.extend(RePar2Less4[i]+np.ones(len(RePar2Less4[i]))*index4)
    index4 += RePar2Less4[i][-1]+1

IsAligned_tot = np.ones(len(RePar2_tot))
IsAlignedLess4_tot = np.ones(len(RePar2Less4_tot))
False_lines = np.zeros(len(P1))

start = -1
for i,chain in zip(range(len(P1Less4)), P1Less4.keys()):
    False_lines[i] = len(P1Less4[chain])+start
    start = False_lines[i]

False_lines = False_lines[:-1]

OverlapandSelfintersectParallelV3(P1Less4_tot, P2Less4_tot, RePar1Less4_tot, RePar2Less4_tot, IsAlignedLess4_tot, P1org, P2org, NresAverage, options, False_lines, P1, P2, RePar1, RePar2, IsAligned,Insert_points_P1, Insert_points_P)





#_______________________________________________________________________________________________________________________





import numpy as np
from NEAMReparametrizationParallel import NEAMReparametrizationParallel
from ScoreSelfIntcWeightedMatchingReparametrizisedParallelTMP import ScoreSelfIntcWeightedMatchingReparametrizisedParallelTMP
from AlignmentMetaData import AlignmentMetaData
from SelfintersectionTransversal import SelfintersectionTransversal
from MakeSelfIntcFigureV3 import MakeSelfIntcFigureV3

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

    Insert_points_P1_tot = np.concatenate(list(Insert_points_P1.values()), axis = 0)
    Insert_points_P_tot = np.concatenate(list(Insert_points_P.values()), axis = 0)
    IPP1_tjek = Insert_points_P1_tot[tjekliste[:, 0]]
    IPP_tjek = Insert_points_P_tot[tjekliste[:, 0]]
    
    tjekliste[:, 0] = tjekliste[:, 0]-IPP1_tjek


    for k in range(tjekliste.shape[0]):
        i = tjekliste[k, 0]
        j = tjekliste[k, 1]
        UdSelf = SelfintersectionTransversal(P1Less4[i:(i+2), :].T, P2Less4[i:(i+2), :].T, P1Less4[j:(j+2), :].T, P2Less4[j:(j+2), :].T)
        UdSelf = np.atleast_2d(UdSelf)
        selfintc[i, j] = UdSelf[0, 0]
        print(f"{k/PotSelfIntc*100:.2f}%")
        if UdSelf[0, 0] ** 2 == 1:
            selfintcu[i, j] = UdSelf[0, 1]
            selfintcv[i, j] = UdSelf[0, 2]
            selfintcs[i, j] = UdSelf[0, 3]
    print(len(np.where(selfintc)[0]))
    
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
            for chain1, chain2 in zip(P1Less4, P2Less4):
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




#_______________________________________________________________________________________________________________________


from Bio.PDB import PDBParser
import Bio.PDB
from Bio.SeqUtils import IUPACData
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB.Polypeptide import PPBuilder, CaPPBuilder
from Bio import Align
from PDBP_to_seq import two_PDB_to_seq, one_PDB_to_seq
from Align_3D import Align_3D
import itertools
import copy

def structural_alignment(pdb_file1, pdb_file2, makefigure = 0):
    
    def find_missing_numbers(arr, n):
        # Calculate the sum of integers from 1 to n
        total_sum = n * (n + 1) // 2
        
        # Calculate the sum of elements in the array
        arr_sum = sum(arr)
        
        # Calculate the difference to find the sum of missing numbers
        # missing_sum = total_sum - arr_sum
        
        # Find the missing numbers
        missing_numbers = []
        for i in range(1, n + 1):
            if i not in arr:
                missing_numbers.append(i)
        
        return missing_numbers
    
    def find_increasing_subarrays(arr):
        # Initialize the current length and the result list
        current_length = 1
        result = []
        result2 = []

        # Iterate over the array
        for i in range(1, len(arr)):
            # If the current number is one greater than the previous number, increase the current length
            if arr[i] == arr[i - 1] + 1:
                current_length += 1
            else:
                # Otherwise, add the current length to the result list current_length times, and reset it
                result.extend(np.linspace(1, current_length, current_length, dtype=int))
                result2.extend([current_length]*current_length)
                current_length = 1

        # Don't forget to add the last subarray length
        result.extend(np.linspace(1, current_length, current_length, dtype=int))
        result2.extend([current_length]*current_length)

        return result, result2
    
    # Function to calculate distance matrix for a permutation
    def distance_matrix_for_permutation(perm):
        best_perm = None
        min_RMSD = np.inf
        best_perms = []
        for letter in perm:
            com_array2 = np.zeros((len(chain_com2), 3))
            i = 0
            for chain in letter:
                for j in range(len(chain_com2[chain])):
                    com_array2[i,j] = chain_com2[chain][j]
                i+=1
            transformed_pts, R, RMSD = Align_3D(np.array(com_array2), np.array(com_array))
            if min_RMSD >= RMSD:
                min_RMSD = RMSD
                best_perms.append(letter)
        best_perms = best_perms[-1]
        return best_perms

    P1, P2, seq1, seq2, ref_structure, sample_structure, tot_seq1, tot_seq2, chain_com1, chain_com2 = two_PDB_to_seq(pdb_file1, pdb_file2)
    P1_org = copy.deepcopy(P1)
    P2_org = copy.deepcopy(P2)


    chain_name1 = list(seq1.keys())
    chain_name2 = list(seq2.keys())


    distance_matrix1 = np.zeros((len(chain_name1), len(chain_name1)))
    distance_matrix2 = np.zeros((len(chain_name2), len(chain_name2)))
    nr_chains = len(chain_name1)

    for i in range(nr_chains):
        for j in range(nr_chains):
            if chain_name1[i] == chain_name1[j]:
                distance_matrix1[i, j] = 0
            else:
                distance_matrix1[i, j] = np.linalg.norm(np.array(chain_com1[chain_name1[i]]) - np.array(chain_com1[chain_name1[j]]))
                distance_matrix1[j, i] = distance_matrix1[i, j]
            if chain_name2[i] == chain_name2[j]:
                distance_matrix2[i, j] = 0
            else:
                distance_matrix2[i, j] = np.linalg.norm(np.array(chain_com2[chain_name2[i]]) - np.array(chain_com2[chain_name2[j]]))
                distance_matrix2[j, i] = distance_matrix2[i, j]

    permutations = list(itertools.permutations(chain_name2))


    i = 0
    com_array = np.zeros((len(chain_com1), 3))

    for chain in chain_com1.keys():
        # Populating the array with values from lists
        for j in range(len(chain_com1[chain])):
            com_array[i,j] = chain_com1[chain][j]
        i += 1

    best_perms = distance_matrix_for_permutation(permutations)



    # Find optimal chain pairs
    Best_chain_pairs = [best_perms]

    #Index for best chain pair
    Best_chain_index = 0

    #Reorder chains in P2 and seq2
    P2_Reorder = {Best_chain_pairs[Best_chain_index][i]: P2[Best_chain_pairs[0][i]] for i in range(len(P2))}
    seq2_Reorder = {Best_chain_pairs[Best_chain_index][i]: seq2[Best_chain_pairs[0][i]] for i in range(len(seq2))}

    chain_name1 = list(seq1.keys())
    chain_name2 = list(seq2_Reorder.keys())



    # Start alignment
    aligner = Align.PairwiseAligner()

    align = {}
    for chain1, chain2 in zip(chain_name1, chain_name2):
        alignments = aligner.align(seq1[chain1], seq2[chain2])
        align[chain1] = alignments[0]
        # print("Score = %.1f:" % alignments[0].score)


    atoms_to_be_aligned1 = {}
    atoms_to_be_aligned2 = {}
    for chain1, chain2 in zip(chain_name1, chain_name2):
        Num_holes = align[chain1].aligned[0].shape[0]
        atoms_to_be_aligned1[chain1] = []
        atoms_to_be_aligned2[chain2] = []
        for i in range(Num_holes-1):
            atoms_to_be_aligned1[chain1].extend(range((align[chain1].aligned[0][i][0]),(align[chain1].aligned[0][i][1])))
            atoms_to_be_aligned2[chain2].extend(range((align[chain1].aligned[1][i][0]),(align[chain1].aligned[1][i][1])))

        atoms_to_be_aligned1[chain1].extend(range((align[chain1].aligned[0][Num_holes-1][0]),(align[chain1].aligned[0][Num_holes-1][1])+1))
        atoms_to_be_aligned2[chain2].extend(range((align[chain1].aligned[1][Num_holes-1][0]),(align[chain1].aligned[1][Num_holes-1][1])+1))



    for chain in P1:
        P1[chain] = P1[chain].tolist()
        P2_Reorder[chain] = P2_Reorder[chain].tolist()

        # Extracting the list of lists from P1
        lists1 = P1[chain]
        lists2 = P2_Reorder[chain]

        # Creating a NumPy array with the same length as lists and 3 columns
        P1_array = np.zeros((len(lists1), 3))
        P2_array = np.zeros((len(lists2), 3))
    
        # Populating the array with values from lists
        for i, sublist in enumerate(lists1):
            P1_array[i] = sublist

        for i, sublist in enumerate(lists2):
            P2_array[i] = sublist

        # Replacing the list of lists with the NumPy array
        P1[chain] = P1_array
        P2_Reorder[chain] = P2_array



    #Center the points
    for chain in P1:
        P1[chain] = P1[chain] - np.mean(P1[chain], axis = 0)
        P2_Reorder[chain] = P2_Reorder[chain] - np.mean(P2_Reorder[chain], axis = 0)

    aligment_points1 = np.zeros((0,3))
    aligment_points2 = np.zeros((0,3))

    for chain1, chain2 in zip(P1, P2_Reorder):
        for i in atoms_to_be_aligned1[chain1]:
            aligment_points1 = np.vstack((aligment_points1, P1[chain1][i-1]))
        for i in atoms_to_be_aligned2[chain2]:
            aligment_points2 = np.vstack((aligment_points2, P2_Reorder[chain2][i-1]))

    aligment_points1 = aligment_points1[1:,:]
    aligment_points2 = aligment_points2[1:,:]


    Transformed_points, R, rmsd = Align_3D(aligment_points1, aligment_points2)

    P = {}
    start = 0
    for chain1, chain2 in zip(P1, P2_Reorder):
        P[chain1] = Transformed_points[start:start+len(atoms_to_be_aligned2[chain1])-1]
        start += len(atoms_to_be_aligned2[chain1])
        # atoms_not_aligned = find_missing_numbers(atoms_to_be_aligned2[chain1], len(P[chain1]))
        #atoms_not_aligned = [i for i, x in enumerate(align[chain1][0]) if x == "-"]
        
        # Find the difference between the two sets
        atoms_not_aligned = set(range(0,len(P1[chain1]))) - set(atoms_to_be_aligned2[chain1])

        # Convert the set to a list
        atoms_not_aligned = sorted(list(atoms_not_aligned))

        for i,j in enumerate(reversed(atoms_not_aligned)):
            P[chain1] = np.insert(P[chain1], j-(5-i), R@P2_Reorder[chain2][j-1], axis=0)

    for chain in P1:
        P1[chain] = P1[chain].tolist()
        P[chain] = P[chain].tolist()

    repar = {}
    repar1 = {}


    for chain in chain_name1:
        repar[chain] = np.linspace(0,len(P[chain])-1,len(P[chain])).tolist()
        repar1[chain] = np.linspace(0,len(P1[chain])-1,len(P1[chain])).tolist()


    indices_target = {}
    indices_query = {}
    
    for key in P:
        indices_target[key] = [i for i, x in enumerate(align[key][1]) if x == "-"]
        indices_query[key]  = [i for i, x in enumerate(align[key][0]) if x == "-"]
        # print(indices_target[key])
        # print(indices_query[key])
        Factor_hole_target, Index_hole_target  = find_increasing_subarrays(indices_target[key])
        Factor_hole_query, Index_hole_query = find_increasing_subarrays(indices_query[key])


        for i in reversed(range(len(indices_target[key]))):
            index = indices_target[key][i]
            alpha = Factor_hole_target[i]/(Index_hole_target[i]+1)
            new_point = [alpha*P[key][index][0]+(1-alpha)*P[key][index+1][0],
                        alpha*P[key][index][1]+(1-alpha)*P[key][index+1][1],
                        alpha*P[key][index][2]+(1-alpha)*P[key][index+1][2]]
            P[key].insert(index+1,new_point)
            repar[key].insert(index+1-(Factor_hole_target[i]-1),index+alpha-(Factor_hole_target[i]-1))
            # print(repar1[key][(index+i-5):(index+i+5)])

        for i in reversed(range(len(indices_query[key]))):
            index = indices_query[key][i]
            alpha = Factor_hole_query[i]/(Index_hole_query[i]+1)
            new_point = [alpha*P1[key][index][0]+(1-alpha)*P1[key][index+1][0],
                        alpha*P1[key][index][1]+(1-alpha)*P1[key][index+1][1],
                        alpha*P1[key][index][2]+(1-alpha)*P1[key][index+1][2]]
            P1[key].insert(index+1,new_point)
            repar1[key].insert(index+1-(Factor_hole_query[i]-1),index+alpha-(Factor_hole_query[i]-1))

    L1 = {}
    L2 = {}
    Insert_points_P1 = {}
    Insert_points_P = {}
    PLess4 = copy.deepcopy(P)
    P1Less4 = copy.deepcopy(P1)

    
    ReParLess4 = copy.deepcopy(repar)
    RePar1Less4 = copy.deepcopy(repar1)
    print("Length of repar[Chain_A]: ", len(repar["Chain_A"]))
    print("Length of repar1[Chain_A]: ", len(repar1["Chain_A"]))
    # Insert points in linesegments  > 4
    for chain1, chain2 in zip(P1Less4, PLess4):
        n = len(P1Less4[chain1])
        m =  len(PLess4[chain2])
        P1_tmp = np.array(P1Less4[chain1])
        P_tmp = np.array(PLess4[chain2])
        L1[chain1] = np.sqrt(np.sum((P1_tmp[0:n - 1, :] - P1_tmp[1:n, :]) ** 2, axis=1))
        L2[chain2] = np.sqrt(np.sum((P_tmp[0:m - 1, :] - P_tmp[1:m, :]) ** 2, axis=1))
        Lmax = np.maximum((L1[chain1]), (L2[chain2]))
        Long_lines = np.where(Lmax > 4)
        Insert_points_P1[chain1] = np.zeros((n)).tolist()
        Insert_points_P[chain2] = np.zeros((m)).tolist()
        
        # print("1: " + str(Long_lines1))
        # print("2: " + str(Long_lines2))
        for i in reversed(Long_lines[0]):
            #print(np.sqrt(np.sum((np.array(P1[chain1])[i,:]-np.array(P1[chain1])[i+1,:])**2)))
            P1Less4[chain1].insert(i+1, ((np.array(P1Less4[chain1])[i,:]+np.array(P1Less4[chain1])[i+1,:])/2).tolist())
            Insert_points_P1[chain1].insert(i+1, 1)
            RePar1Less4[chain1].insert(i+1, (RePar1Less4[chain1][i]+RePar1Less4[chain1][i+1])/2)

            # print(max(np.sqrt(np.sum((P_tmp[i,:]-P_tmp[i+1,:])**2)),np.sqrt(np.sum((P1_tmp[i,:]-P1_tmp[i+1,:])**2))))
            PLess4[chain2].insert(i+1, ((np.array(PLess4[chain2])[i,:]+np.array(PLess4[chain2])[i+1,:])/2).tolist())
            Insert_points_P[chain2].insert(i+1, 1)
            ReParLess4[chain2].insert(i+1, (ReParLess4[chain2][i]+ReParLess4[chain2][i+1])/2)

    print("Length of repar[Chain_A]: ", len(repar["Chain_A"]))
    print("Length of repar1[Chain_A]: ", len(repar1["Chain_A"]))



    # Lav repar
    if makefigure == 1:
        # #Plot P1, P2 and P in 3d using plotly
        import plotly.graph_objects as go

        fig = go.Figure()

        for chain in P1.keys():
            fig.add_trace(go.Scatter3d(x=[i[0] for i in P1[chain]], y=[i[1] for i in P1[chain]], z=[i[2] for i in P1[chain]], mode='lines', line=dict(width=9, color = "blue"), name=chain))

        # for chain in P2.keys():
        #     fig.add_trace(go.Scatter3d(x=[i[0] for i in P2[chain]], y=[i[1] for i in P2[chain]], z=[i[2] for i in P2[chain]], mode='lines', line=dict(width=9), name=chain))

        for chain in P.keys():
            fig.add_trace(go.Scatter3d(x=[i[0] for i in P[chain]], y=[i[1] for i in P[chain]], z=[i[2] for i in P[chain]], mode='lines', line=dict(width=9,color = 'red'), name="Aligned "+chain))

        #add plot title
        fig.update_layout(title_text="Structural alignment of protein structures")
        fig.show()

        pv1 = 270
        pv2 = 285

        
        #Create a plot for each pair of chains
        for i in range(len(P1.keys())):
            fig = go.Figure()
            fig.add_trace(go.Scatter3d(x=[i[0] for i in P[chain_name1[i]][pv1:pv2]], y=[i[1] for i in P[chain_name1[i]][pv1:pv2]], z=[i[2] for i in P[chain_name1[i]][pv1:pv2]], mode='lines', line=dict(width=9), name='P'))
            fig.add_trace(go.Scatter3d(x=[i[0] for i in P[chain_name1[i]][pv1:pv2]], y=[i[1] for i in P[chain_name1[i]][pv1:pv2]], z=[i[2] for i in P[chain_name1[i]][pv1:pv2]], mode='lines', line=dict(width=9), name='P'))
            fig.update_layout(title_text="Structural alignment of protein structures for chain " + chain_name1[i])
            fig.show()

    print("RMSD of structual alignment " + str(rmsd))
    # print(best_perms)

    is_aligned = {}
    NresAverage = {}

    for chain in repar:
        is_aligned[chain] = np.ones(len(repar1[chain]))
        # NresAverage[chain] = (len(P1_org[chain])+len(P2_org[chain]))/2
        P1[chain] = np.array(P1[chain])
        P[chain] = np.array(P[chain])

    P1org_tot = np.concatenate(list(P1_org.values()), axis = 0)
    P2org_tot = np.concatenate(list(P2_org.values()), axis = 0)
    # print(P2_org)
    NresAverage = (len(P1org_tot)+len(P2org_tot))/2


    return P1, P, repar1, repar, is_aligned, NresAverage, P1Less4, PLess4, RePar1Less4, ReParLess4, Insert_points_P1, Insert_points_P


#pdb_file1 = "/Users/agb/Desktop/Bachelor projekt/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRUA_hexamer_positive.pdb"
#pdb_file2 = "/Users/agb/Desktop/Bachelor projekt/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRU1_hexamer_negative.pdb"

#pdb_file1 = "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB//CRUA_hexamer_positive.pdb"
#pdb_file2 = "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRU1_hexamer_negative.pdb"


#P1, P, repar1, repar, is_aligned, NresAverage = structural_alignment(pdb_file1, pdb_file2, makefigure = 1)


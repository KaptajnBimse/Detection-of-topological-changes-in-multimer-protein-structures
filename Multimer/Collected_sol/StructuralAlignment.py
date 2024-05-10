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
from FindIncreasingSubarray import find_increasing_subarrays
from DistanceMatrix import distance_matrix_permutations
import itertools
import copy


def structural_alignment(pdb_file1, pdb_file2, makefigure = 0):
    P1, P2, seq1, seq2, chain_com1, chain_com2 = two_PDB_to_seq(pdb_file1, pdb_file2)
    
    # Create copies of each configuration
    P1_org = copy.deepcopy(P1)
    P2_org = copy.deepcopy(P2)
    
    # Obtain list of chain names
    chain_name1 = list(seq1.keys())
    chain_name2 = list(seq2.keys())

    # Create all possible permutations of chain names
    permutations = list(itertools.permutations(chain_name2))
    com_array = np.zeros((len(chain_com1), 3))
    
    i = 0
    for chain in chain_com1.keys():
        for j in range(len(chain_com1[chain])):
            com_array[i,j] = chain_com1[chain][j]
        i += 1

    # Find best permutation of chain names
    best_perms = distance_matrix_permutations(permutations,chain_com2,com_array)
    best_chain_pairs = [best_perms]
    
    # Choose first permutation (change index number if another of the best permutations is desired)
    best_chain_index = 0

    #Reorder chains in P2 and seq2 to new optimal order
    P2_Reorder = {best_chain_pairs[best_chain_index][i]: P2[best_chain_pairs[0][i]] for i in range(len(P2))}
    seq2_Reorder = {best_chain_pairs[best_chain_index][i]: seq2[best_chain_pairs[0][i]] for i in range(len(seq2))}

    chain_name1 = list(seq1.keys())
    chain_name2 = list(seq2_Reorder.keys())

    # Start alignment
    aligner = Align.PairwiseAligner()

    align = {}
    for chain1, chain2 in zip(chain_name1, chain_name2):
        alignments = aligner.align(seq1[chain1], seq2[chain2])
        align[chain1] = alignments[0]

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

    # Convert coordinate list to single array for each chain
    for chain in P1:
        P1[chain] = P1[chain].tolist()
        P2_Reorder[chain] = P2_Reorder[chain].tolist()

        # Extracting the coordinates of residues from P1
        lists1 = P1[chain]
        lists2 = P2_Reorder[chain]
        P1_array = np.zeros((len(lists1), 3))
        P2_array = np.zeros((len(lists2), 3))
    
        # Convert the lists to arrays
        for i, sublist in enumerate(lists1):
            P1_array[i] = sublist

        for i, sublist in enumerate(lists2):
            P2_array[i] = sublist

        # Store the arrays in the dictionaries
        P1[chain] = P1_array
        P2_Reorder[chain] = P2_array

    #Center coordinates of each chain around 0
    for chain in P1:
        P1[chain] = P1[chain] - np.mean(P1[chain], axis = 0)
        P2_Reorder[chain] = P2_Reorder[chain] - np.mean(P2_Reorder[chain], axis = 0)

    alignment_points1 = np.zeros((0,3))
    alignment_points2 = np.zeros((0,3))

    # Create connected total structure
    for chain1, chain2 in zip(P1, P2_Reorder):
        for i in atoms_to_be_aligned1[chain1]:
            alignment_points1 = np.vstack((alignment_points1, P1[chain1][i-1]))
        for i in atoms_to_be_aligned2[chain2]:
            alignment_points2 = np.vstack((alignment_points2, P2_Reorder[chain2][i-1]))

    alignment_points1 = alignment_points1[1:,:]
    alignment_points2 = alignment_points2[1:,:]

    # Align the two connected structures
    Transformed_points, R, rmsd = Align_3D(alignment_points1, alignment_points2)

    P = {}
    start = 0
    for chain1, chain2 in zip(P1, P2_Reorder):
        P[chain1] = Transformed_points[start:start+len(atoms_to_be_aligned2[chain1])-1]
        start += len(atoms_to_be_aligned2[chain1])
        
        # Find the difference between the two sets
        atoms_not_aligned = set(range(0,len(P1[chain1]))) - set(atoms_to_be_aligned2[chain1])

        # Convert the set to a list
        atoms_not_aligned = sorted(list(atoms_not_aligned))

        for i,j in enumerate(reversed(atoms_not_aligned)):
            P[chain1] = np.insert(P[chain1], j-(5-i), R@P2_Reorder[chain2][j-1], axis=0)

    # Create RePar
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
    
    # Insert interpolation points in the aligned structures
    for key in P:
        indices_target[key] = [i for i, x in enumerate(align[key][1]) if x == "-"]
        indices_query[key]  = [i for i, x in enumerate(align[key][0]) if x == "-"]

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
    
    # Initiate copies for inserting points in linesegments > 4
    PLess4 = copy.deepcopy(P)
    P1Less4 = copy.deepcopy(P1)
    ReParLess4 = copy.deepcopy(repar)
    RePar1Less4 = copy.deepcopy(repar1)

    
    # Insert points in linesegments > 4
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
        
        for i in reversed(Long_lines[0]):
            P1Less4[chain1].insert(i+1, ((np.array(P1Less4[chain1])[i,:]+np.array(P1Less4[chain1])[i+1,:])/2).tolist())
            Insert_points_P1[chain1].insert(i+1, 1)
            RePar1Less4[chain1].insert(i+1, (RePar1Less4[chain1][i]+RePar1Less4[chain1][i+1])/2)

            PLess4[chain2].insert(i+1, ((np.array(PLess4[chain2])[i,:]+np.array(PLess4[chain2])[i+1,:])/2).tolist())
            Insert_points_P[chain2].insert(i+1, 1)
            ReParLess4[chain2].insert(i+1, (ReParLess4[chain2][i]+ReParLess4[chain2][i+1])/2)

    if makefigure == 1:
        # #Plot P1, P2 and P in 3d using plotly
        import plotly.graph_objects as go

        fig = go.Figure()

        for chain in P1.keys():
            fig.add_trace(go.Scatter3d(x=[i[0] for i in P1[chain]], y=[i[1] for i in P1[chain]], z=[i[2] for i in P1[chain]], mode='lines', line=dict(width=9, color = "blue"), name=chain))

        for chain in P.keys():
            fig.add_trace(go.Scatter3d(x=[i[0] for i in P[chain]], y=[i[1] for i in P[chain]], z=[i[2] for i in P[chain]], mode='lines', line=dict(width=9,color = 'red'), name="Aligned "+chain))
            
        fig.update_layout(title_text="Structural alignment of protein structures")
        fig.show()

        pv1 = 0

        #Create a plot for each pair of chains
        for i in range(len(P1.keys())):
            pv2 = len(P1[chain_name1[i]])
            fig = go.Figure()
            fig.add_trace(go.Scatter3d(x=[i[0] for i in P[chain_name1[i]][pv1:pv2]], 
                                       y=[i[1] for i in P[chain_name1[i]][pv1:pv2]], 
                                       z=[i[2] for i in P[chain_name1[i]][pv1:pv2]], mode='lines', line=dict(width=9), name='P'))
            
            fig.add_trace(go.Scatter3d(x=[i[0] for i in P[chain_name1[i]][pv1:pv2]], y=[i[1] for i in P[chain_name1[i]][pv1:pv2]], z=[i[2] for i in P[chain_name1[i]][pv1:pv2]], mode='lines', line=dict(width=9), name='P'))
            fig.update_layout(title_text="Structural alignment of protein structures for chain " + chain_name1[i])
            fig.show()

    print("RMSD of structual alignment " + str(rmsd))

    is_aligned = {}
    NresAverage = {}

    for chain in repar:
        is_aligned[chain] = np.ones(len(repar1[chain]))
        P1[chain] = np.array(P1[chain])
        P[chain] = np.array(P[chain])

    # Create org totals
    P1org_tot = np.concatenate(list(P1_org.values()), axis = 0)
    P2org_tot = np.concatenate(list(P2_org.values()), axis = 0)
    NresAverage = (len(P1org_tot)+len(P2org_tot))/2

    return P1, P, repar1, repar, is_aligned, NresAverage, P1Less4, PLess4, RePar1Less4, ReParLess4, Insert_points_P1, Insert_points_P
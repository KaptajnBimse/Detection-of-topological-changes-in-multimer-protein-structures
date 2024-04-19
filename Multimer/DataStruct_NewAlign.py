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


#P1, P2, seq1, seq2, ref_structure, sample_structure, tot_seq1, tot_seq2= two_PDB_to_seq("examples/Multimer PDB/CRUA_hexamer_positive.pdb", "examples/Multimer PDB/CRU1_hexamer_negative.pdb")
P1, P2, seq1, seq2, ref_structure, sample_structure, tot_seq1, tot_seq2, chain_com1, chain_com2= two_PDB_to_seq("C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRUA_hexamer_positive.pdb"
, "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRU1_hexamer_negative.pdb"
)

# Find optimal chain pairs
Best_chain_pairs = [('Chain_A', 'Chain_B', 'Chain_C', 'Chain_D', 'Chain_E', 'Chain_F'), ('Chain_D', 'Chain_C', 'Chain_B', 'Chain_A', 'Chain_E', 'Chain_F'), ('Chain_B', 'Chain_A', 'Chain_D', 'Chain_C', 'Chain_E', 'Chain_F'), ('Chain_B', 'Chain_A', 'Chain_D', 'Chain_C', 'Chain_F', 'Chain_E')]
#Index for best chain pair
Best_chain_index = 1
#Reorder chains in P2 and seq2
P2_Reorder = {Best_chain_pairs[Best_chain_index][i]: P2[Best_chain_pairs[0][i]] for i in range(len(P2))}
seq2_Reorder = {Best_chain_pairs[Best_chain_index][i]: seq2[Best_chain_pairs[0][i]] for i in range(len(seq2))}

chain_name1 = list(seq1.keys())
chain_name2 = list(seq2.keys())


aligner = Align.PairwiseAligner()

align = {}
for chain in seq1:
    alignments = aligner.align(seq1[chain], seq2[chain])
    align[chain] = alignments[0]
    print("Score = %.1f:" % alignments[0].score)


atoms_to_be_aligned = {}
for chain in seq1:
  Num_holes = align[chain].aligned[0].shape[0]
  atoms_to_be_aligned[chain] = []
  for i in range(Num_holes):
    atoms_to_be_aligned[chain].extend(range((align[chain].aligned[0][i][0]+1),(align[chain].aligned[0][i][1])+1))



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

for chain in P1:
    for i in atoms_to_be_aligned[chain]:
        aligment_points1 = np.vstack((aligment_points1, P1[chain][i-1]))
        aligment_points2 = np.vstack((aligment_points2, P2_Reorder[chain][i-1]))

aligment_points1 = aligment_points1[1:,:]
aligment_points2 = aligment_points2[1:,:]


Transformed_points, R, rmsd = Align_3D(aligment_points1, aligment_points2)

P = {}
start = 0
for chain in P1:
    P[chain] = Transformed_points[start:start+len(atoms_to_be_aligned[chain])]
    start += len(atoms_to_be_aligned[chain])

for chain in P1:
    P1[chain] = P1[chain].tolist()
    P[chain] = P[chain].tolist()

#Kode der fjerener eventuel ekstra punkter start eller slut i structure
#...
#-------

repar = {}
repar1 = {}


for chain in chain_name2:
    repar[chain] = np.linspace(0,len(P[chain]),len(P[chain])+1).tolist()
    repar1[chain] = np.linspace(0,len(P1[chain]),len(P1[chain])+1).tolist()


indices_target = {}
indices_query = {}
for key in P:
  indices_target[key] = [i for i, x in enumerate(align[key][0]) if x == "-"]
  indices_query[key]  = [i for i, x in enumerate(align[key][1]) if x == "-"]

  Factor_hole_target, Index_hole_target  = find_increasing_subarrays(indices_target[key])
  Factor_hole_query, Index_hole_query = find_increasing_subarrays(indices_query[key])

  for i in range(len(indices_target[key])):
    index = indices_target[key][i]
    alpha = Index_hole_target[i]/(Factor_hole_target[i]+1)
    new_point = [alpha*P1[key][index-(Index_hole_target[i]-1)][0]+(1-alpha)*P1[key][index-(Index_hole_target[i]-1)+1][0],
                alpha*P1[key][index-(Index_hole_target[i]-1)][1]+(1-alpha)*P1[key][index-(Index_hole_target[i]-1)+1][1],
                alpha*P1[key][index-(Index_hole_target[i]-1)][2]+(1-alpha)*P1[key][index-(Index_hole_target[i]-1)+1][2]]
    P1[key].insert(index+1,new_point)
    repar1[key].insert(index+i,index+alpha)
    print(repar1[key][(index+i-5):(index+i+5)])

  for i in range(len(indices_query[key])):
    index = indices_query[key][i]
    alpha = Index_hole_query[i]/(Factor_hole_query[i]+1)
    new_point = [alpha*P[key][index-(Index_hole_query[i]-1)][0]+(1-alpha)*P[key][index-(Index_hole_query[i]-1)+1][0],
                alpha*P[key][index-(Index_hole_query[i]-1)][1]+(1-alpha)*P[key][index-(Index_hole_query[i]-1)+1][1],
                alpha*P[key][index-(Index_hole_query[i]-1)][2]+(1-alpha)*P[key][index-(Index_hole_query[i]-1)+1][2]]
    P[key].insert(index+1,new_point)
    repar[key].insert(index+1,index+alpha)

# Lav repar



# #Plot P1, P2 and P in 3d using plotly
import plotly.graph_objects as go


# Truep1 = np.loadtxt("Multimer/Test txt/TestAlign/P1.txt")
# Truep2 = np.loadtxt("Multimer/Test txt/TestAlign/P2.txt")

# trace1 = go.Scatter3d(
#         x=Truep1[:, 0],
#         y=Truep1[:, 1],
#         z=Truep1[:, 2],
#         mode='lines',
#         line=dict(width=9),
#         name='True 1',
#         legendgroup= 'Chain 1',
#     )

# trace2 = go.Scatter3d(
#         x=Truep2[:, 0],
#         y=Truep2[:, 1],
#         z=Truep2[:, 2],
#         mode='lines',
#         line=dict( width=9),
#         name='True 1',
#         legendgroup= 'Chain 1',
#     )

fig = go.Figure()

for chain in P1.keys():
    fig.add_trace(go.Scatter3d(x=[i[0] for i in P1[chain]], y=[i[1] for i in P1[chain]], z=[i[2] for i in P1[chain]], mode='lines', line=dict(width=9, color = "red"), name=chain))

# for chain in P2.keys():
#     fig.add_trace(go.Scatter3d(x=[i[0] for i in P2[chain]], y=[i[1] for i in P2[chain]], z=[i[2] for i in P2[chain]], mode='lines', line=dict(width=9), name=chain))

for chain in P.keys():
    fig.add_trace(go.Scatter3d(x=[i[0] for i in P[chain]], y=[i[1] for i in P[chain]], z=[i[2] for i in P[chain]], mode='lines', line=dict(width=9,color = 'blue'), name="Aligned "+chain))




fig.show()

#Create a plot for each pair of chains

for i in range(len(P1.keys())):
    fig = go.Figure()
    fig.add_trace(go.Scatter3d(x=[i[0] for i in P1[chain_name1[i]]], y=[i[1] for i in P1[chain_name1[i]]], z=[i[2] for i in P1[chain_name1[i]]], mode='lines', line=dict(width=9), name='P1'))
    fig.add_trace(go.Scatter3d(x=[i[0] for i in P[chain_name2[i]]], y=[i[1] for i in P[chain_name2[i]]], z=[i[2] for i in P[chain_name2[i]]], mode='lines', line=dict(width=9), name='P'))
    fig.show()

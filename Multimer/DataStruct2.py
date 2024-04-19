from Bio.PDB import PDBParser
import Bio.PDB
from Bio.SeqUtils import IUPACData
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB.Polypeptide import PPBuilder, CaPPBuilder
from Bio import Align
from PDBP_to_seq import two_PDB_to_seq, one_PDB_to_seq
import itertools

def find_increasing_subarrays(arr):
    # Initialize the current length and max length
    current_length = 1
    lengths = []

    # Iterate over the array
    for i in range(1, len(arr)):
        # If the current number is one greater than the previous number, increase the current length
        if arr[i] == arr[i - 1] + 1:
            current_length += 1
        else:
            # Otherwise, store the current length and reset it
            lengths.append(current_length)
            current_length = 1

    # Don't forget to add the last subarray length
    lengths.append(current_length)

    return lengths


#P1, P2, seq1, seq2, ref_structure, sample_structure, tot_seq1, tot_seq2= two_PDB_to_seq("examples/Multimer PDB/CRUA_hexamer_positive.pdb", "examples/Multimer PDB/CRU1_hexamer_negative.pdb")
P1, P2, seq1, seq2, ref_structure, sample_structure, tot_seq1, tot_seq2, chain_com1, chain_com2 = two_PDB_to_seq("/Users/agb/Desktop/Bachelor projekt/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRUA_hexamer_positive.pdb"
, "/Users/agb/Desktop/Bachelor projekt/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRU1_hexamer_negative.pdb"
)

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

# row1 = np.sum(distance_matrix1, axis = 0)
# row2 = np.sum(distance_matrix2, axis = 0)
# min_array = np.full((len(row1), 3), np.inf)
# for i in range(len(row1)):
#     for j in range(len(row2)):
#         if j not in min_array[:,2]:
#             min_temp = np.abs(row1[i]-row2[j])
#             print(min_temp)
#             if min_temp < min_array[i,0]:
#                 min_array[i,0] = min_temp
#                 min_array[i,1] = i
#                 min_array[i,2] = j

# permutations = list(itertools.permutations(chain_name2))

# # Function to calculate distance matrix for a permutation
# def distance_matrix_for_permutation(perm):
#     best_perm = None
#     sz = np.inf
#     for letter in perm:
#         for i in range(nr_chains):
#             for j in range(nr_chains):
#                 if chain_name2[i] == chain_name2[j]:
#                     distance_matrix2[i, j] = 0
#                 else:
#                     distance_matrix2[i, j] = np.linalg.norm(np.array(chain_com2[letter[i]]) - np.array(chain_com2[letter[j]]))
#                     distance_matrix2[j, i] = distance_matrix2[i, j]
#         diff = np.abs(distance_matrix1 - distance_matrix2)
#         if sz > np.linalg.norm(diff):
#             sz = np.linalg.norm(diff)
#             best_perm = letter
#     return best_perm
# best_perm = distance_matrix_for_permutation(permutations)

permutations = list(itertools.permutations(chain_name2))

# Function to calculate distance matrix for a permutation
def distance_matrix_for_permutation(perm,tol):
    best_perm = None
    sz = np.inf
    best_perms = []
    for letter in perm:
        for i in range(nr_chains):
            for j in range(nr_chains):
                if chain_name2[i] == chain_name2[j]:
                    distance_matrix2[i, j] = 0
                else:
                    distance_matrix2[i, j] = np.linalg.norm(np.array(chain_com2[letter[i]]) - np.array(chain_com2[letter[j]]))
                    distance_matrix2[j, i] = distance_matrix2[i, j]
        diff = np.abs(distance_matrix1 - distance_matrix2)
        norm_diff = np.linalg.norm(diff)
        if (1/(1+tol))*sz > norm_diff:
            sz = norm_diff
            best_perms.append(letter)
    return best_perms

best_perms = distance_matrix_for_permutation(permutations,0.1)


# bestemme massemidtpunkt for hver kæde i query og target hexamer.
# TARGET = den vi ønsker at fitte til den originale
# danne afstandsmatrice som euklidisk afstand mellem massemidtpunkt for kæder i.
# fastholde target afstandsmatrice

# lave permutationer af matricen og tjekke euklidisk afstand mellem to afstandsmatricer
# hvis denne afstand ligger inden for en vis afstand til hinanden (brugerdefineret) skal vi beholde flere af permutationerne (vi kan godt rotere hexameren i flere klik som hver især passer rimeligt godt allesammen).
# tag det antal permutationsmatricer med videre til alignment og superimposer og se hvad de giver
# hvis RMSD er inden for bestemt interval her også - tag med videre igen
# tjek til sidst hvilken morph der så er bedst


aligner = Align.PairwiseAligner()

align = {}
# hvis flere sekvenser lav for-loop udenom nedenstående for-loop
# seq1 skal ændres til den optimale
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


# Select what residues numbers you wish to align
# and put them in a list
# start_id = align.aligned[1][0][0]+1
# end_id   = align.aligned[1][0][1]
# atoms_to_be_aligned = range(start_id, end_id + 1)

# Use the first model in the pdb-files for alignment
# Change the number 0 if you want to align to another structure
ref_model    = ref_structure[0]
sample_model = sample_structure[0]

# Make a list of the atoms (in the structures) you wish to align.
# In this case we use CA atoms whose index is in the specified range
ref_atoms = []
sample_atoms = []

# Iterate of all chains in the model in order to find all residues
index = 0
for ref_chain in ref_model:
  # Iterate of all residues in each model in order to find proper atoms
  for ref_res in ref_chain:
    # Check if residue number ( .get_id() ) is in the list
    if ref_res.get_id()[1] in atoms_to_be_aligned[chain_name1[index]]:
      # Append CA atom to list
      # print(ref_res.get_id()[1])
      ref_atoms.append(ref_res['CA'])
  index += 1


# Do the same for the sample structure
index = 0
for sample_chain in sample_model:
  for sample_res in sample_chain:
    if sample_res.get_id()[1] in atoms_to_be_aligned[chain_name2[index]]:
      # print(sample_res.get_id()[1])
      sample_atoms.append(sample_res['CA'])
  index += 1

# Now we initiate the superimposer:
super_imposer = Bio.PDB.Superimposer()
super_imposer.set_atoms(ref_atoms, sample_atoms)
super_imposer.apply(sample_model.get_atoms())

# Print RMSD:
print(super_imposer.rms)

# Save the aligned version of 1UBQ.pdb
io = Bio.PDB.PDBIO()
io.set_structure(sample_structure) 
aligned = io.structure


ca = [(atom.full_id[2], atom.full_id[3][1], *atom.get_coord()) for atom in aligned.get_atoms() if atom.full_id[4][0] == "CA"]
df = pd.DataFrame(ca, columns=["chain", "residue_number", "x", "y", "z"])


# x_mean = df["x"].mean()
# y_mean = df["y"].mean()
# z_mean = df["z"].mean()

# df["x"] = df["x"] - x_mean
# df["y"] = df["y"] - y_mean
# df["z"] = df["z"] - z_mean

df["P"] = df[["x", "y", "z"]].values.tolist()

#Kode der fjerener eventuel ekstra punkter start eller slut i structure
#...
#-------

P = {}
for chain in df["chain"].unique():
    P["Chain_" + chain] = df[df["chain"] == chain]["P"].array


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
    fig.add_trace(go.Scatter3d(x=[i[0] for i in P1[chain]], y=[i[1] for i in P1[chain]], z=[i[2] for i in P1[chain]], mode='lines', line=dict(width=9), name=chain))

# for chain in P2.keys():
#     fig.add_trace(go.Scatter3d(x=[i[0] for i in P2[chain]], y=[i[1] for i in P2[chain]], z=[i[2] for i in P2[chain]], mode='lines', line=dict(width=9), name=chain))

for chain in P.keys():
    fig.add_trace(go.Scatter3d(x=[i[0] for i in P[chain]], y=[i[1] for i in P[chain]], z=[i[2] for i in P[chain]], mode='lines', line=dict(width=9), name="Aligned "+chain))

fig.show()

#Create a plot for each pair of chains
for chain in P1.keys():
    fig = go.Figure()
    fig.add_trace(go.Scatter3d(x=[i[0] for i in P1[chain]], y=[i[1] for i in P1[chain]], z=[i[2] for i in P1[chain]], mode='lines', line=dict(width=9), name='P1'))
    fig.add_trace(go.Scatter3d(x=[i[0] for i in P[chain]], y=[i[1] for i in P[chain]], z=[i[2] for i in P[chain]], mode='lines', line=dict(width=9), name='P'))
    fig.show()

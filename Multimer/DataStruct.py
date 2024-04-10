from Bio.PDB import PDBParser
import Bio.PDB
from Bio.SeqUtils import IUPACData
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB.Polypeptide import PPBuilder
from Bio import Align
from PDBP_to_seq import two_PDBP_to_seq

P1, P2, seq1, seq2, ref_structure, sample_structure= two_PDBP_to_seq("C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/T1026.pdb", "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/T1026TS015_5")

aligner = Align.PairwiseAligner()

alignments = aligner.align(seq1["Chain_ "], seq2["Chain_ "])
# print(len(alignments))
# for alignment in sorted(alignments):

#     print("Score = %.1f:" % alignment.score)

#     print(alignment)

align = sorted(alignments)[0]

# Select what residues numbers you wish to align
# and put them in a list
start_id = align.aligned[1][0][0]+1
end_id   = align.aligned[1][0][1]
atoms_to_be_aligned = range(start_id, end_id + 1)

# Use the first model in the pdb-files for alignment
# Change the number 0 if you want to align to another structure
ref_model    = ref_structure[0]
sample_model = sample_structure[0]

# Make a list of the atoms (in the structures) you wish to align.
# In this case we use CA atoms whose index is in the specified range
ref_atoms = []
sample_atoms = []

# Iterate of all chains in the model in order to find all residues
for ref_chain in ref_model:
  # Iterate of all residues in each model in order to find proper atoms
  for ref_res in ref_chain:
    # Check if residue number ( .get_id() ) is in the list
    if ref_res.get_id()[1] in atoms_to_be_aligned:
      # Append CA atom to list
      ref_atoms.append(ref_res['CA'])


# Do the same for the sample structure
for sample_chain in sample_model:
  for sample_res in sample_chain:
    if sample_res.get_id()[1] in atoms_to_be_aligned:
      sample_atoms.append(sample_res['CA'])

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


x_mean = df["x"].mean()
y_mean = df["y"].mean()
z_mean = df["z"].mean()

df["x"] = df["x"] - x_mean
df["y"] = df["y"] - y_mean
df["z"] = df["z"] - z_mean

df["P"] = df[["x", "y", "z"]].values.tolist()


P = {}
for chain in df["chain"].unique():
    P["Chain_" + chain] = df[df["chain"] == chain]["P"]



#Plot P1, P2 and P in 3d using plotly
import plotly.graph_objects as go


Truep1 = np.loadtxt("Multimer/Test txt/TestAlign/P1.txt")
Truep2 = np.loadtxt("Multimer/Test txt/TestAlign/P2.txt")

trace1 = go.Scatter3d(
        x=Truep1[:, 0],
        y=Truep1[:, 1],
        z=Truep1[:, 2],
        mode='lines',
        line=dict(width=9),
        name='True 1',
        legendgroup= 'Chain 1',
    )

trace2 = go.Scatter3d(
        x=Truep2[:, 0],
        y=Truep2[:, 1],
        z=Truep2[:, 2],
        mode='lines',
        line=dict( width=9),
        name='True 1',
        legendgroup= 'Chain 1',
    )

fig = go.Figure(data=[trace2])

for chain in P1.keys():
    fig.add_trace(go.Scatter3d(x=[i[0] for i in P1[chain]], y=[i[1] for i in P1[chain]], z=[i[2] for i in P1[chain]], mode='lines', line=dict(width=9), name=chain))

# for chain in P2.keys():
#     fig.add_trace(go.Scatter3d(x=[i[0] for i in P2[chain]], y=[i[1] for i in P2[chain]], z=[i[2] for i in P2[chain]], mode='lines', line=dict(width=9), name=chain))

# for chain in P.keys():
#     fig.add_trace(go.Scatter3d(x=[i[0] for i in P[chain]], y=[i[1] for i in P[chain]], z=[i[2] for i in P[chain]], mode='lines', line=dict(width=9), name="Aligned "+chain))




fig.show()
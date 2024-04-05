from Bio.PDB import PDBParser
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

p = PDBParser(QUIET=True)
s = p.get_structure("s", "examples/Multimer PDB/1y8k.pdb")

ca = [(atom.full_id[2], atom.full_id[3][1], *atom.get_coord()) for atom in s.get_atoms() if atom.full_id[4][0] == "CA"]
df = pd.DataFrame(ca, columns=["chain", "residue_number", "x", "y", "z"])
#df

df["P"] = df[["x", "y", "z"]].values.tolist()

P = {}

for chain in df["chain"].unique():
    P["Chain_" + chain] = df[df["chain"] == chain]["P"]


from Bio.pairwise2 import nwalign

# Define two sequences
seq1 = "ACCGT"
seq2 = "ACG"

# Perform Needleman-Wunsch alignment
alignments = nwalign(seq1, seq2)

# Print the alignment results
for alignment in alignments:
    print("Alignment Score:", alignment.score)
    print("Sequence 1:", alignment.seqA)
    print("Sequence 2:", alignment.seqB)
    print("Alignment Start:", alignment.start)
    print("Alignment End:", alignment.end)
    print()
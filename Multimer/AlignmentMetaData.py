import numpy as np

def AlignmentMetaData(RePar1, RePar2, IsAligned):
    """
    AlignmentMetaData returns [#aligned residues, length of aligned window 
    chain 1, length of aligned window chain 2, # vertices in the full alignment]
    """
    out = [np.sum(IsAligned), RePar1[-1] - RePar1[0] + 1,
           RePar2[-1] - RePar2[0] + 1, len(RePar1)]
    return out

# Repar1 = np.loadtxt('Test txt/AlignmentMetaData/RePar1.txt')
# Repar2 = np.loadtxt('Test txt/AlignmentMetaData/RePar2.txt')
# IsAligned = np.loadtxt('Test txt/AlignmentMetaData/IsAligned.txt')
# print(AlignmentMetaData(Repar1, Repar2, IsAligned))

""" 
from Bio import pairwise2
from Bio.PDB import PDBParser

# Parse PDB files
parser = PDBParser()
structure1 = parser.get_structure('pdb1', 'pdb1_filename.pdb')
structure2 = parser.get_structure('pdb2', 'pdb2_filename.pdb')

# Extract sequences
chain1 = structure1[0]['A']
chain2 = structure2[0]['A']
seq1 = ''.join([residue.get_resname() for residue in chain1])
seq2 = ''.join([residue.get_resname() for residue in chain2])

# Perform alignment
alignments = pairwise2.align.globalxx(seq1, seq2)

# Print the alignment
for alignment in alignments:
    print(pairwise2.format_alignment(*alignment))
    
    

# Hexamer:
from Bio import pairwise2
from Bio.PDB import PDBParser

# Parse PDB files
parser = PDBParser()
structure1 = parser.get_structure('pdb1', 'pdb1_filename.pdb')
structure2 = parser.get_structure('pdb2', 'pdb2_filename.pdb')

# Extract sequences for all chains in each hexamer
hexamer1_seq = ''
for chain in structure1[0]:
    hexamer1_seq += ''.join([residue.get_resname() for residue in chain])

hexamer2_seq = ''
for chain in structure2[0]:
    hexamer2_seq += ''.join([residue.get_resname() for residue in chain])

# Perform alignment
alignments = pairwise2.align.globalxx(hexamer1_seq, hexamer2_seq)

# Print the alignment
for alignment in alignments:
    print(pairwise2.format_alignment(*alignment))

# This code assumes that each hexamer is represented as the first model ([0]) in each structure. 




# Using StructureAlignment
from Bio.PDB import PDBParser, StructureAlignment

# Parse PDB files
parser = PDBParser()
structure1 = parser.get_structure('pdb1', 'pdb1_filename.pdb')
structure2 = parser.get_structure('pdb2', 'pdb2_filename.pdb')

# Create StructureAlignment object
alignment = StructureAlignment()
alignment.add_structure('pdb1', structure1)
alignment.add_structure('pdb2', structure2)

# Print RMSD and aligned residues
print(alignment)
# This will print the RMSD (Root Mean Square Deviation) and aligned residues between the two structures. 
# You can access the aligned residues and RMSD values for further analysis.





# Using superimposer
from Bio.PDB import PDBParser, Superimposer

# Parse PDB files
parser = PDBParser()
structure1 = parser.get_structure('pdb1', 'pdb1_filename.pdb')
structure2 = parser.get_structure('pdb2', 'pdb2_filename.pdb')

# Superimpose the structures
superimposer = Superimposer()
atoms1 = []
atoms2 = []
for atom1, atom2 in zip(structure1.get_atoms(), structure2.get_atoms()):
    atoms1.append(atom1)
    atoms2.append(atom2)
superimposer.set_atoms(atoms1, atoms2)

# Print RMSD and transformation matrix
print("RMSD:", superimposer.rms)
print("Rotation matrix:\n", superimposer.rotran[0])
print("Translation vector:", superimposer.rotran[1])
# This will print the RMSD (Root Mean Square Deviation) 
# between the two structures and the transformation matrix needed to superimpose one structure onto the other. """
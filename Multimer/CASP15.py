import numpy as np
from Structural_AlignmentV2 import structural_alignment
from OverlapandSelfintersectParallelV3 import OverlapandSelfintersectParallelV3
from PDBP_to_seq import two_PDB_to_seq, one_PDB_to_seq
from FinalFunction import FinalFunction
import os



# pdb_file1 = "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CASP15T1132/T1132TS180_1o.pdb"
# pdb_file2 = "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CASP15T1132/T1132TS462_1o.pdb"
# pdb_file1 = "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRUA_hexamer_positive.pdb"
# pdb_file2 = "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRU1_hexamer_negative.pdb"
#one_PDB_to_seq("C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CASP15Target/T1124o.pdb")
target = "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CASP15Target/T1109o.pdb"

# Specify the directory you want to scan
directory = "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/T1109o/"

# Get a list of all files in the directory
files = os.listdir(directory)

# Filter the list to include only .pdb files
pdb_files = [f for f in files]

# Run the function on each .pdb file
for pdb_file in pdb_files:
    full_path = os.path.join(directory, pdb_file)
    FinalFunction(target, full_path)
    print("Finished", pdb_file)
# FinalFunction(pdb_file1, pdb_file2)

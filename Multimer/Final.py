import numpy as np
from Structural_Alignment import structural_alignment
from OverlapandSelfintersectParallelV3 import OverlapandSelfintersectParallelV3

pdb_file1 = "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB//CRUA_hexamer_positive.pdb"
pdb_file2 = "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRU1_hexamer_negative.pdb"


P1, P2, RePar1, RePar2, IsAligned, NresAverage =  structural_alignment(pdb_file1, pdb_file2)
options = {'Smoothning': 1, 'AllowEndContractions': 1, 'MaxLength': 1, 'MakeFigures': 1}
P1org = 0
P2org = 0

OverlapandSelfintersectParallelV3(P1, P2, RePar1, RePar2, IsAligned, P1org, P2org, NresAverage, options)
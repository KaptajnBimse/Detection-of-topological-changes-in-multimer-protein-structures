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

OverlapandSelfintersectParallelV3(P1Less4_tot, P2Less4_tot, RePar1Less4_tot, RePar2Less4_tot, IsAlignedLess4_tot, P1org, P2org, NresAverage, options, False_lines, P1, P2, RePar1_tot, RePar2_tot, IsAligned,Insert_points_P1, Insert_points_P)

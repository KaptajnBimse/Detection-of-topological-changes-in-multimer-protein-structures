import numpy as np
from Structural_AlignmentV2 import structural_alignment
from OverlapandSelfintersectParallelV3 import OverlapandSelfintersectParallelV3

Adam = 0
if Adam == 1:
    pdb_file1 = "/Users/agb/Desktop/Bachelorprojekt/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRUA_hexamer_positive.pdb"
    pdb_file2 = "/Users/agb/Desktop/Bachelorprojekt/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRU1_hexamer_negative.pdb"
else:
    pdb_file1 = "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRUA_hexamer_positive.pdb"
    pdb_file2 = "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRU1_hexamer_negative.pdb"

    # pdb_file1 = "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/1y8h.pdb"
    # pdb_file2 = "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/1y8i.pdb"
    

P1, P2, RePar1, RePar2, IsAligned, NresAverage, P1Less4, P2Less4, RePar1Less4, RePar2Less4, Insert_points_P1, Insert_points_P, b_factors1, b_factors2 =  structural_alignment(pdb_file1, pdb_file2, makefigure = 0)
# options = {'Smoothning': 0, 'AllowEndContractions': 0, 'MaxLength': 5, 'MakeFigures': 1}
options = {
    'MaxLength': 10,
    'dmax': 10,
    'Smoothning': 0,
    'AllowEndContractions': 0,
    'MakeFigures': 1,
    'MakeAlignmentSeedFigure': 0,
    'MakeFiguresInLastItteration': 1,
    'MakeLocalPlotsOfEssensials': 1,
    'SelfIntcFigCutSize': 10,
    'PrintOut': 0,
    'additionalRMSD': 0,
    'alignmentsmoothing': 0,
    'alignmentsmoothingwidth': 3,
    'AdaptiveSubset': 1,
    'MaxNbrAlignmentSeeds': 7,
    'MaxSeedOverlap': 0.5000,
    'MinSeedLength': 40,
    'OverlapWeight': 4,
    'MaxIter': 20,
    'MaxWindowMisalignment': 1,
    'MaxMisAlignment': 0.0150,
    'MinimalAlignmentLength': 30,
    'FileName1': 'CRUA_hexamer_positive.pdb',
    'FileName2': 'CRUA_hexamer_negative.pdb',
    'StructureSequenceWeight': 1.5608,
    'SeqenceMisAlignmentPenalty': [7.2200  ,  2.1660], 
    'TrimSeqenceAlignment': 0,
    'SequenceAlignmentExtension': 1,
    'InitialAlignmentExactPairs': 1
}
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

OverlapandSelfintersectParallelV3(P1Less4_tot, P2Less4_tot, RePar1Less4_tot, RePar2Less4_tot, IsAlignedLess4_tot, P1org, P2org, NresAverage, options, False_lines, P1, P2, RePar1_tot, RePar2_tot, IsAligned,Insert_points_P1, Insert_points_P, b_factors1, b_factors2)

"""
'calls': Sort by call count.
'cumulative': Sort by cumulative time.
'filename': Sort by the name of the file in which the function was defined.
'line': Sort by the line number in the file where the function was defined.
'module': Sort by the name of the module in which the function was defined.
'name': Sort by function name.
'nfl': Sort by name, file, and line number.
'pcalls': Sort by primitive call count.
'stdname': Sort by standard name.
'time': Sort by internal time.
"""

import cProfile
import pstats
import matplotlib.pyplot as plt
import os
from collections import defaultdict

# def profile_and_print_stats(func, *args, **kwargs):
#     profiler = cProfile.Profile()
#     profiler.enable()
#     func(*args, **kwargs)
#     profiler.disable()
    
#     stats = pstats.Stats(profiler).sort_stats('time')
    
#     script_times = defaultdict(float)
#     for func_name, info in stats.stats.items():
#         filename = func_name[0]
#         # replace 'your_script_names' with the names of your scripts
#         if any(script in filename for script in ['OverlapandSelfintersectParallelV3.py', 'AlignmentMetaData.py', 'NEAMReparametrizationParallel','SelfintersectionTransversal','ScoreSelfIntcWeightedMatchingReparametrizisedParallelTMP','MakeSelfIntcFigureV3' 'Final.py']):
#             script_times[os.path.basename(filename)] += info[2]  # total time
    
#     plt.barh(list(script_times.keys()), list(script_times.values()), color='blue')
#     plt.xlabel('Total Time')
#     plt.ylabel('Script')
#     plt.title('Script Execution Time')
#     plt.show()

# def my_function(P1Less4_tot, P2Less4_tot, RePar1Less4_tot, RePar2Less4_tot, IsAlignedLess4_tot, P1org, P2org, NresAverage, options, False_lines, P1, P2, RePar1_tot, RePar2_tot, IsAligned,Insert_points_P1, Insert_points_P):
#     OverlapandSelfintersectParallelV3(P1Less4_tot, P2Less4_tot, RePar1Less4_tot, RePar2Less4_tot, IsAlignedLess4_tot, P1org, P2org, NresAverage, options, False_lines, P1, P2, RePar1_tot, RePar2_tot, IsAligned,Insert_points_P1, Insert_points_P)

# profile_and_print_stats(my_function,P1Less4_tot, P2Less4_tot, RePar1Less4_tot, RePar2Less4_tot, IsAlignedLess4_tot, P1org, P2org, NresAverage, options, False_lines, P1, P2, RePar1_tot, RePar2_tot, IsAligned,Insert_points_P1, Insert_points_P)







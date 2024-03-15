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
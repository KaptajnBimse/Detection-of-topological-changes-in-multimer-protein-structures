import numpy as np
from Align_3D import Align_3D

# Function to calculate distance matrix for a permutation
def distance_matrix_permutations(perm,chain_com2,com_array):
    best_perm = None
    min_RMSD = np.inf
    best_perms = []
    for letter in perm:
        com_array2 = np.zeros((len(chain_com2), 3))
        i = 0
        for chain in letter:
            for j in range(len(chain_com2[chain])):
                com_array2[i,j] = chain_com2[chain][j]
            i+=1
        transformed_pts, R, RMSD = Align_3D(np.array(com_array2), np.array(com_array))
        if min_RMSD >= RMSD:
            min_RMSD = RMSD
            best_perms.append(letter)
    best_perms = best_perms[-1]
    return best_perms

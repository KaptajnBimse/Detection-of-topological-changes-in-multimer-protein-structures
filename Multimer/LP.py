# import numpy as np
# from scipy.optimize import linear_sum_assignment

# # Example: distances_matrix_1 and distances_matrix_2 are distance matrices between center of mass points for two structures
# distances_matrix_1 = np.array([[0,3,1,2], [3,0,5,1], [1,5,0,2], [2,1,2,0]])  # n1 and n2 are the number of chains in the two structures
# distances_matrix_2 = np.array([[0,a,b,c], [d,0,e,f], [g,h,0,i], [j,k,l,0]])

# # Calculate cost matrix for linear sum assignment (scipy.optimize.linear_sum_assignment)
# cost_matrix = np.abs(distances_matrix_1 - distances_matrix_2)

# # Solve linear sum assignment problem to find best alignment
# row_indices, col_indices = linear_sum_assignment(cost_matrix)

# # Align chains based on the assignment
# aligned_chains = [(row_idx, col_idx) for row_idx, col_idx in zip(row_indices, col_indices)]

# # Output aligned chains
# print("Aligned chains:")
# for row_idx, col_idx in aligned_chains:
#     print(f"Chain {row_idx} in Structure 1 aligns with Chain {col_idx} in Structure 2")


import numpy as np
from itertools import permutations
# lave permutationer af bogstaverne (kædenavn) og lave distance matrix ud fra det
def permutation_matrices(n):
    identity_matrix = np.eye(n)
    row_permutations = permutations(identity_matrix)
    permutation_matrices = [np.array(perm) for perm in row_permutations]
    return permutation_matrices

# Generate permutation matrices for a 6x6 matrix
permutation_matrices_6x6 = permutation_matrices(4)

# Print the permutation matrices
for idx, perm_matrix in enumerate(permutation_matrices_6x6, start=1):
    print(f"Permutation Matrix {idx}:")
    print(perm_matrix)
    print()


import itertools
points = chain_com2
letters = ['A', 'B', 'C', 'D', 'E', 'F']
permutations = list(itertools.permutations(letters))
def distance_matrix_for_permutation(perm):
    permutation_points = [points["Chain_"+letter] for letter in perm]
    distances = np.linalg.norm(permutation_points - np.expand_dims(permutation_points, axis=1), axis=2)
    return distances

import numpy as np
import numpy.matlib

def MakeDP(P):
    # P is assumed to be an n times 3 matrix
    n, three = P.shape # n is the number of points, three is number of coordinates
    if three != 3: # If the number of coordinates is not 3, return an empty array
        dP = np.array([])
        return dP
    
    dP = np.zeros((n, n, 3)) # Create an n times n times 3 array
    for i in range(3):
        dP[:, :, i] = np.tile(P[:, i], (n, 1)) - np.tile(P[:, i], (n, 1)).T
    return dP

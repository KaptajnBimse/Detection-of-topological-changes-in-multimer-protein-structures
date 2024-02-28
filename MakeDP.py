import numpy as np

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

# Example usage
#np.random.seed(123)
#PP = np.random.rand(10, 3)  # Generate a random matrix P with 1000 rows and 3 columns
#dP = MakeDP(PP)  # Call MakeDP function with input matrix P
#print(dP)  # Print the resulting dP array

# P1 = np.array([
#     [-15.9862, 0.5931, -26.8763],
#     [-18.0812, -0.3186, -23.7314],
#     [-20.4852, -3.2195, -24.5687],
#     [-22.4635, -2.2352, -21.6092],
#     [-20.1591, -2.7425, -18.6242],
#     [-20.8995, -1.7523, -14.9922],
#     [-17.9512, 0.2156, -13.6301],
#     [-19.7483, 0.9920, -10.3647],
#     [-18.9081, -2.5654, -9.0747],
#     [-17.9987, -4.3071, -12.3942],
#     [-16.3278, -7.5281, -11.5771],
#     [-15.1077, -9.9705, -14.3713],
#     [-15.9844, -7.3638, -16.8812],
#     [-13.0136, -6.8066, -19.0663]
# ])
# MakeDP(P1)
# P2 = np.array([
#     [-13.0512, 4.3497, -24.4231],
#     [-10.7082, 1.7397, -22.9551],
#     [-12.6772, -1.0183, -21.1551],
#     [-12.2912, -0.1983, -17.4631],
#     [-12.9652, -2.5443, -14.5611],
#     [-12.8592, -1.6903, -10.8611],
#     [-11.8382, -4.8413, -9.0301],
#     [-11.4752, -5.2183, -5.2721],
#     [-8.6622, -7.1613, -3.6361],
#     [-8.1962, -7.4713, 0.1229],
#     [-5.0202, -8.3893, 1.9829],
#     [-5.0042, -9.2743, 5.6689],
#     [-1.6932, -8.8913, 7.4689],
#     [-0.9662, -9.1373, 11.1739]
# ])
# DP = MakeDP(P2)
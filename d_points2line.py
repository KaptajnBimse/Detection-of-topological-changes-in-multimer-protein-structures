import numpy as np

def d_points2line(Ps, P1, P2):
    V = P2 - P1
    V = V / np.sqrt(np.sum(V ** 2))
    n = Ps.shape[1]
    Vs = np.tile(V.reshape(-1, 1), (1, n))
    vcross = np.cross(Ps - np.tile(P1.reshape(-1, 1), (1, n)), Vs, axis=0)
    ud = np.sqrt(np.sum(vcross ** 2, axis=0))
    return ud

# P1 = np.loadtxt("C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Python code/Oversæt/Test txt/d_points2line/P1.txt")
# P2 = np.loadtxt("C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Python code/Oversæt/Test txt/d_points2line/P2.txt")
# Ps = np.loadtxt("C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Python code/Oversæt/Test txt/d_points2line/Ps.txt")

# d_points2line(Ps, P1, P2)
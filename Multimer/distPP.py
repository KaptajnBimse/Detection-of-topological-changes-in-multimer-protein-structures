import numpy as np

def distPP(p1, p2):
    ud = np.sqrt(np.sum((p2 - p1) ** 2))
    return ud
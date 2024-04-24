import numpy as np
import MakeDminProteinReparametrizedParallel as mdprp
import MakeDP as mdp


def NEAMReparametrizationParallel(P1, P2, RePar1, RePar2, IsAligned, Smoothning):

    Dmin = mdprp.MakeDminProteinReparametrizedParallel(RePar1, RePar2, Smoothning)
    dP1 = mdp.MakeDP(P1)
    dP2 = mdp.MakeDP(P2)
    l1 = np.sqrt(np.sum(dP1**2, axis=2))
    ddPNormSqr = np.sum((dP1 - dP2)**2, axis=2)
    dot12 = np.sum(dP1 * dP2, axis=2)
    crossPNormSqr = np.sum(np.cross(dP2, dP1, axis=2)**2, axis=2)
    dminSqr = crossPNormSqr / ddPNormSqr
    t = (l1**2 - dot12) / ddPNormSqr
    tstar = np.maximum(np.minimum(t, 1), 0)
    sEffSq = (tstar - t)**2
    dminSqrSegment = sEffSq * ddPNormSqr + dminSqr
    overlap = Dmin - np.sqrt(dminSqrSegment)
    overlapalt = Dmin - np.sqrt(np.sum(dP1**2, axis=2))
    overlap[ddPNormSqr < 1.0E-15] = overlapalt[ddPNormSqr < 1.0E-15]
    overlap = np.maximum(overlap, 0)
<<<<<<< HEAD
    #___________________________________________________________________________
=======

    
>>>>>>> 9e03e534733566c26700b52e47caa4253520f42c
    alignedaligned = np.outer(IsAligned, IsAligned)
    overlapaligned = overlap * alignedaligned
    tmp1 = np.diff(RePar1)
    weight1 = 0.5 * (np.concatenate(([1], tmp1)) + np.concatenate((tmp1, [1])))
    tmp2 = np.diff(RePar2)
    weight2 = 0.5 * (np.concatenate(([1], tmp2)) + np.concatenate((tmp2, [1])))
    overlapGap = overlap - overlapaligned
    overlapGapWeight = overlapGap * np.outer(weight1, weight2)

    return overlap, overlapaligned, overlapGap, overlapGapWeight

# P1 = data = np.loadtxt("C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Python code/Oversæt/P1.txt")
# P2 = data = np.loadtxt("C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Python code/Oversæt/P2.txt")
# RePar1 = data = np.loadtxt("C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Python code/Oversæt/RePar1.txt")
# RePar2 = data = np.loadtxt("C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Python code/Oversæt/RePar2.txt")
# IsAligned = data = np.loadtxt("C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Python code/Oversæt/IsAligned.txt")

# overlap, overlapaligned, overlapGap, overlapGapWeight = NEAMReparametrizationParallel(P1, P2, RePar1, RePar2, IsAligned, 1)


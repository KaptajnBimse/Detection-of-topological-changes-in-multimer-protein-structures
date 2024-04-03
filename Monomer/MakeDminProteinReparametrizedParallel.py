import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import splrep, PPoly

def MakeDminProteinReparametrizedParallel(RePar1, RePar2, Smoothning):
    # Returns the minimal distance allowed between residues i,j as function
    # of the smaller of the "arclength distance" between the i'th and j'th
    # point in the two reparametrizations.

    SCALEFACTOR = 1  # 1
    nbr_points = len(RePar1)
    if abs(nbr_points - len(RePar2)):
        print('Error: MakeDminProteinReparametrized requiers equal length parametrizations')
        return

    if Smoothning == 1:
        tmp = np.array([1.0, 2.1, 3.0, 3.4, 3.6, 3.7, 3.7])
        mind = 3.7 * np.ones(nbr_points)
        mind[:min(7, nbr_points)] = tmp[:min(7, nbr_points)]
    else:
        tmp = np.array([2.8, 4.5, 3.86, 3.47, 3.52, 3.48, 3.6])
        mind = 3.7 * np.ones(nbr_points)
        mind[:min(7, nbr_points)] = tmp[:min(7, nbr_points)]
    mind = mind * SCALEFACTOR
    # sp = BSpline(np.arange(nbr_points + 2), mind, k=3)
    # tck = splrep(np.arange(nbr_points), mind)
    tck = splrep(np.arange(nbr_points+1), np.hstack((0,mind)), k = 3)
    pp = PPoly.from_spline(tck)

    tmp = np.tile(RePar1, (nbr_points)).reshape(nbr_points, nbr_points)
    diff_arclength1 = np.abs(tmp - tmp.T)
    tmp2 = np.tile(RePar2, (nbr_points)).reshape(nbr_points, nbr_points)
    diff_arclength2 = np.abs(tmp2 - tmp2.T)
    # Dminalt = sp(np.minimum(diff_arclength1, diff_arclength2))
    Dminalt = pp(np.minimum(diff_arclength1, diff_arclength2))
    return Dminalt

# RePar1 = np.arange(1, 17) + 23

# RePar2 = np.array([1, 1.2, 1.4, 1.6, 1.8, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
# Dminalt,da1,da2 = MakeDminProteinReparametrizedParallel(RePar1,RePar2,1)
# print(Dminalt[0:3,0:3])
# import matplotlib.pyplot as plt
# plt.plot(np.minimum(da1,da2)[0,:],Dminalt[0,:])
# plt.show()

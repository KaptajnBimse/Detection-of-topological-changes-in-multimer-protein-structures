import numpy as np
from scipy.linalg import lstsq

def intersection_triangle_line_segment(p0, p1, p2, Lstart, Lslut):
    A = np.array([p1 - p0, p2 - p0, -Lslut + Lstart]).T
    b = (Lstart - p0)
    uvt,_,_,_ = np.linalg.lstsq(A, b, rcond=None)
    uvt = uvt.reshape(-1,1)
    ud = uvt[0,0] >= 0 and uvt[1,0] >= 0 and uvt[0,0] + uvt[1,0] <= 1 and uvt[2,0] >= 0 and uvt[2,0] <= 1
    ud = int(ud)
    return ud, uvt

# p0 = np.array([10.5149, 0.9082, -4.0328])
# p1 = np.array([12.3474, -0.9595, -4.5598])
# p2 = np.array([10.9492, 1.1307, -3.7315])
# Lstart = np.array([12.0179,-0.8118,-5.2360])
# Lslut = np.array([11.9545,-1.2143,-5.0977])

# [ud, uvt] = intersection_triangle_line_segment(p0, p1, p2, Lstart, Lslut)


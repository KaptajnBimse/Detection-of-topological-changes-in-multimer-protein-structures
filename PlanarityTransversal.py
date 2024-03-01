import numpy as np

def PlanarityTransversal(a0, a1, b0, b1):
    # a(t) = (1-t)*[a_0^0 a_0^1] + t*[a_1^0 a_1^1]
    # b(t) = (1-t)*[b_0^0 b_0^1] + t*[b_1^0 b_1^1]
    # These are functions of t that move the endpoints of the linesegments 
    # (i.e. a0 is the line segment at t=0 and a1 is the line segment at t=1)
    # The function returns the values of t where the linesegments intersect (that is, the a and b lines intersect)
    # The line segments are here just made from a 3*2 matrix, since a line is constructed between 2 end points
    
    # The push of the line segment can be seen as a vector v_a(t) and v_b(t) that moves the line segments
    # To find where these two line segments are planar, we need to find the values of t where the line segments intersect 
    # (however this is too complicated since we end up with the lines L1(s_1,t) = L2(s_2,t) and this calculation is too much)
    # Insteead we use the following method:
    # We find the determinant of the matrix M = [v_a(t) v_b(t)] and the vector k = b_0(t) - a_0(t)
    # This means when the determinant is zero, the lines are parallel and we can find the intersection by finding the determinant of the cross product of the vectors
    # It will be a polynomial of degree 3 in t, with real roots t in ]0,1[
    
    # a0 is 3x2 and holds the a-linesegment at s=0
    # a1 is 3x2 and holds the a-linesegment at s=1
    # similarly for b0 and b1
    a = a0[:,0]  # point on first line segment at s=0
    b = a0[:,1] - a  # vector in first linesegment at s=0
    da = a1[:,0] - a  # change of a form s=0 to s=1
    db = a1[:,1] - a1[:,0] - b  # do for b
    c = b0[:,0]  # now for the second linesegment
    d = b0[:,1] - c
    dc = b1[:,0] - c
    dd = b1[:,1] - b1[:,0] - d
    if (np.dot(da, da) + np.dot(db, db) + np.dot(dc, dc) + np.dot(dd, dd)) < 1.0e-25:
        ud = [[], []]  # both linsegments do not move
        return ud
    r = c - a
    dr = dc - da
    M1 = np.column_stack((dr, r, dr, r))
    M2 = np.column_stack((db, db, b, b))
    d_matrix = np.array([np.transpose(d),np.transpose(dd)])
    M_matrix = np.array([M1[1, :] * M2[2, :] - M1[2, :] * M2[1, :],
                         M1[2, :] * M2[0, :] - M1[0, :] * M2[2, :],
                         M1[0, :] * M2[1, :] - M1[1, :] * M2[0, :]])
    M = np.dot(d_matrix, M_matrix)

    p = np.array([M[1, 0], M[0, 0] + M[1, 1] + M[1, 2], M[0, 1] + M[0, 2] + M[1, 3], M[0, 3]])
    tmp = np.roots(p)
    dp = [3, 2, 1] * p[:3]
    tmp = tmp[np.imag(tmp) == 0]
    tmp = tmp[tmp >= 0]
    ud = [tmp[tmp <= 1], []] if tmp.size > 0 else [[], []]
    if len(ud[0]) > 0:
        ud[1] = np.dot(np.column_stack((ud[0] ** 2, ud[0], np.ones(len(ud[0])))), dp)
    return ud

a0 = np.loadtxt("/Users/agb/Desktop/Bachelor projekt/Python kode oversat/a0.txt")
a1 = np.loadtxt("/Users/agb/Desktop/Bachelor projekt/Python kode oversat/a1.txt")
b0 = np.loadtxt("/Users/agb/Desktop/Bachelor projekt/Python kode oversat/b0.txt")
b1 = np.loadtxt("/Users/agb/Desktop/Bachelor projekt/Python kode oversat/b1.txt")

udplan = PlanarityTransversal(a0,a1,b0,b1)



import numpy as np
import PlanarityTransversal as PT
import plotly.graph_objects as go

def SelfintersectionTransversal(a0, a1, b0, b1):
    ud = 0
    udplan = PT.PlanarityTransversal(a0, a1, b0, b1)
    slist = udplan[0] # t values of planarity
    transversal = udplan[1] # derivative of planes at these t values
    cut = 10**-20
    i = 0
    #for s in slist:
    for s in (x for x in slist if x > 0):
        # i+= 1
        a = (1 - s) * a0 + s * a1
        b = (1 - s) * b0 + s * b1
        M = np.column_stack((a[:, 1] - a[:, 0], b[:, 0] - b[:, 1]))
        k = b[:, 0] - a[:, 0]
        if np.linalg.matrix_rank(M) == 1:
            tmp = np.sum(np.cross(M[:, 0], k)**2) # The length of the vector between the parallel lines (orthogonal vector to first line)
            if tmp > 10**-20:
                return ud
            else: # Scenario where the lines lay on top of each other
                v = M[:, 0]
                ta = np.dot(v, a)
                tb = np.dot(v, b)
                intersectionlength = abs(ta[0] - tb[1]) + abs(ta[1] - tb[0]) - abs(ta[0] - tb[0]) - abs(ta[1] - tb[1])
                if intersectionlength > 4 * 10**-14:
                    ud = [np.sign(transversal[i, 0]), [0.5, 0.5], s] # Fjern muligvis [] inde i []
                    return ud
                else:
                    return ud
        uv, _, _, _ = np.linalg.lstsq(M, k, rcond=None)
        if np.sum((0 <= uv) & (uv <= 1)) == 2:
            if np.sum((M @ uv - k)**2) < cut:
                ud = [np.sign(transversal[i]), uv[0],uv[1], s] # ændret fra transversal[i,0] fordi andet index var ulovligt, men skal måske bruges

                # fig = go.Figure()
                # fig.add_trace(go.Scatter3d(x=((1-s)*a0+s*a1)[0,:].tolist(), 
                #                         y=((1-s)*a0+s*a1)[1,:].tolist(), 
                #                         z=((1-s)*a0+s*a1)[2,:].tolist(), 
                #                         mode='lines', line=dict(width=9,color = 'red'), name='line 1'))
                
                # fig.add_trace(go.Scatter3d(x=((1-s)*b0+s*b1)[0,:].tolist(), 
                #                         y=((1-s)*b0+s*b1)[1,:].tolist(), 
                #                         z=((1-s)*b0+s*b1)[2,:].tolist(), 
                #                         mode='lines', line=dict(width=9,color = 'blue'), name='line 2'))
                
                # fig.show()
            return ud
        i += 1
    return ud
# a0 = np.loadtxt("/Users/agb/Desktop/Bachelor projekt/Python kode oversat/a0.txt")
# a1 = np.loadtxt("/Users/agb/Desktop/Bachelor projekt/Python kode oversat/a1.txt")
# b0 = np.loadtxt("/Users/agb/Desktop/Bachelor projekt/Python kode oversat/b0.txt")
# b1 = np.loadtxt("/Users/agb/Desktop/Bachelor projekt/Python kode oversat/b1.txt")

# ud=SelfintersectionTransversal(a0,a1,b0,b1)
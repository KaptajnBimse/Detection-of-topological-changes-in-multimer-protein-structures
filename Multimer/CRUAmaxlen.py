import numpy as np
from Structural_AlignmentV2 import structural_alignment
from OverlapandSelfintersectParallelV3 import OverlapandSelfintersectParallelV3
from PDBP_to_seq import two_PDB_to_seq, one_PDB_to_seq
from FinalFunction import FinalFunction
import os
import numpy as np
from Structural_AlignmentV2 import structural_alignment
from OverlapandSelfintersectParallelV3 import OverlapandSelfintersectParallelV3
from PDBP_to_seq import two_PDB_to_seq, one_PDB_to_seq
from FinalFunction import FinalFunction
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go



pdb_file1 = "/Users/agb/Desktop/Bachelorprojekt/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRUA_hexamer_positive.pdb"
pdb_file2 = "/Users/agb/Desktop/Bachelorprojekt/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRU1_hexamer_negative.pdb"

max_len = np.arange(1,51,1)

i = 0
Num_ess = np.zeros(len(max_len))

for maxlen in max_len:
    ud = FinalFunction(pdb_file1, pdb_file2, maxlen)
    Num_ess[i] = ud[0].shape[0]
    print("A maximal length of ", maxlen, " gives this number of essential intersections: "+ str(Num_ess[i]))
    i += 1

#plot the number of essential intersections
import matplotlib.pyplot as plt
plt.plot(Num_ess)
plt.show()

df = pd.DataFrame({'Num_ess': Num_ess, 'Max Length': max_len[-1]})

# Create a line plot
fig = go.Figure(data=go.Scatter(x=list(range(len(Num_ess))), y=Num_ess))

# Customize the layout
fig.update_layout(
    xaxis_title="Max Length",
    yaxis_title="Number of essential intersections"
)

fig.show()




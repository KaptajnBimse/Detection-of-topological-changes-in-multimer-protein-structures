import numpy as np
from Structural_AlignmentV2 import structural_alignment
from OverlapandSelfintersectParallelV3 import OverlapandSelfintersectParallelV3
from PDBP_to_seq import two_PDB_to_seq, one_PDB_to_seq
from FinalFunction import FinalFunction
import os
import pandas as pd
import plotly.express as px

# pdb_file1 = "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CASP15T1132/T1132TS180_1o.pdb"
# pdb_file2 = "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CASP15T1132/T1132TS462_1o.pdb"
# pdb_file1 = "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRUA_hexamer_positive.pdb"
# pdb_file2 = "C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CRU1_hexamer_negative.pdb"
#one_PDB_to_seq("C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/CASP15Target/T1124o.pdb")
target = "/Users/agb/Desktop/Bachelorprojekt/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/T1187o.pdb"

# Specify the directory you want to scan
directory = "/Users/agb/Desktop/Bachelorprojekt/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/T1187o/"

# Get a list of all files in the directory
files = os.listdir(directory)

# Only use files that are not .txt files
pdb_files = sorted([f for f in files if f != '.DS_Store'])
print("")
# Run the function on each .pdb file
i = 0
Num_ess = np.zeros(len(pdb_files))
Num_inter = np.zeros(len(pdb_files))
Num_ess2 = np.zeros(len(pdb_files))
Num_inter2 = np.zeros(len(pdb_files))

maxlen = 18
for pdb_file in pdb_files:
    full_path = os.path.join(directory, pdb_file)
    ud = FinalFunction(target, full_path,maxlen)
    Num_ess[i] = ud[0].shape[0]
    Num_inter[i] = ud[1]
    print("Finished " + pdb_file + " Number of essential intersections "+ str(Num_ess[i]))
    i += 1

maxlen = 50
j = 0
for pdb_file in pdb_files:
    full_path = os.path.join(directory, pdb_file)
    ud2 = FinalFunction(target, full_path,maxlen)
    Num_ess2[j] = ud2[0].shape[0]
    Num_inter2[j] = ud2[1]
    print("Finished " + pdb_file + " Number of essential intersections "+ str(Num_ess[j]))
    j += 1



#Load txt file

TM_score = np.loadtxt( "/Users/agb/Downloads/TM_Score_true.txt", delimiter=",")

#Create data frame with number of essential intersections and TM score
df = pd.DataFrame({'Num_ess': Num_ess, 'TM-score': TM_score, 'Num_inter': Num_inter})
df2 = pd.DataFrame({'Num_ess': Num_ess2, 'TM-score': TM_score, 'Num_inter': Num_inter2})


Ess = np.zeros(6)
OK = np.zeros(6)
No_Ess = np.zeros(6)
Ess2 = np.zeros(6)
OK2 = np.zeros(6)
No_Ess2 = np.zeros(6)

OK[0] = len(df[(df['TM-score'] < 0.5) & (df['Num_ess'] == 0)])
OK[1] = len(df[(df['TM-score'] >= 0.5) & (df['TM-score'] < 0.6) & (df['Num_ess'] == 0)])
OK[2] = len(df[(df['TM-score'] >= 0.6) & (df['TM-score'] < 0.7) & (df['Num_ess'] == 0)])
OK[3] = len(df[(df['TM-score'] >= 0.7) & (df['TM-score'] < 0.8) & (df['Num_ess'] == 0)])
OK[4] = len(df[(df['TM-score'] >= 0.8) & (df['TM-score'] < 0.9) & (df['Num_ess'] == 0)])
OK[5] = len(df[(df['TM-score'] >= 0.9) & (df['Num_ess'] == 0)])

OK2[0] = len(df2[(df2['TM-score'] < 0.5) & (df2['Num_ess'] == 0)])
OK2[1] = len(df2[(df2['TM-score'] >= 0.5) & (df2['TM-score'] < 0.6) & (df2['Num_ess'] == 0)])
OK2[2] = len(df2[(df2['TM-score'] >= 0.6) & (df2['TM-score'] < 0.7) & (df2['Num_ess'] == 0)])
OK2[3] = len(df2[(df2['TM-score'] >= 0.7) & (df2['TM-score'] < 0.8) & (df2['Num_ess'] == 0)])
OK2[4] = len(df2[(df2['TM-score'] >= 0.8) & (df2['TM-score'] < 0.9) & (df2['Num_ess'] == 0)])
OK2[5] = len(df2[(df2['TM-score'] >= 0.9) & (df2['Num_ess'] == 0)])

No_Ess[0] = len(df[(df['TM-score'] < 0.5) & (df['Num_ess'] == 0) & (df['Num_inter'] > 10)])
No_Ess[1] = len(df[(df['TM-score'] >= 0.5) & (df['TM-score'] < 0.6) & (df['Num_ess'] == 0) & (df['Num_inter'] > 10)])
No_Ess[2] = len(df[(df['TM-score'] >= 0.6) & (df['TM-score'] < 0.7) & (df['Num_ess'] == 0) & (df['Num_inter'] > 10)])
No_Ess[3] = len(df[(df['TM-score'] >= 0.7) & (df['TM-score'] < 0.8) & (df['Num_ess'] == 0) & (df['Num_inter'] > 10)])
No_Ess[4] = len(df[(df['TM-score'] >= 0.8) & (df['TM-score'] < 0.9) & (df['Num_ess'] == 0) & (df['Num_inter'] > 10)])
No_Ess[5] = len(df[(df['TM-score'] >= 0.9) & (df['Num_ess'] == 0) & (df['Num_inter'] > 10)])

No_Ess2[0] = len(df2[(df2['TM-score'] < 0.5) & (df2['Num_ess'] == 0) & (df2['Num_inter'] > 10)])
No_Ess2[1] = len(df2[(df2['TM-score'] >= 0.5) & (df2['TM-score'] < 0.6) & (df2['Num_ess'] == 0) & (df2['Num_inter'] > 10)])
No_Ess2[2] = len(df2[(df2['TM-score'] >= 0.6) & (df2['TM-score'] < 0.7) & (df2['Num_ess'] == 0) & (df2['Num_inter'] > 10)])
No_Ess2[3] = len(df2[(df2['TM-score'] >= 0.7) & (df2['TM-score'] < 0.8) & (df2['Num_ess'] == 0) & (df2['Num_inter'] > 10)])
No_Ess2[4] = len(df2[(df2['TM-score'] >= 0.8) & (df2['TM-score'] < 0.9) & (df2['Num_ess'] == 0) & (df2['Num_inter'] > 10)])
No_Ess2[5] = len(df2[(df2['TM-score'] >= 0.9) & (df2['Num_ess'] == 0) & (df2['Num_inter'] > 10)])

#Find rows where tm_score is between 0.5 and 0.6 and number of essential intersections is greater than 0
Ess[0] = len(df[(df['TM-score'] < 0.5) & (df['Num_ess'] > 0)])
Ess[1] = len(df[(df['TM-score'] >= 0.5) & (df['TM-score'] < 0.6) & (df['Num_ess'] > 0)])
Ess[2] = len(df[(df['TM-score'] >= 0.6) & (df['TM-score'] < 0.7) & (df['Num_ess'] > 0)])
Ess[3] = len(df[(df['TM-score'] >= 0.7) & (df['TM-score'] < 0.8) & (df['Num_ess'] > 0)])
Ess[4] = len(df[(df['TM-score'] >= 0.8) & (df['TM-score'] < 0.9) & (df['Num_ess'] > 0)])
Ess[5] = len(df[(df['TM-score'] >= 0.9) & (df['Num_ess'] > 0)])

Ess2[0] = len(df2[(df2['TM-score'] < 0.5) & (df2['Num_ess'] > 0)])
Ess2[1] = len(df2[(df2['TM-score'] >= 0.5) & (df2['TM-score'] < 0.6) & (df2['Num_ess'] > 0)])
Ess2[2] = len(df2[(df2['TM-score'] >= 0.6) & (df2['TM-score'] < 0.7) & (df2['Num_ess'] > 0)])
Ess2[3] = len(df2[(df2['TM-score'] >= 0.7) & (df2['TM-score'] < 0.8) & (df2['Num_ess'] > 0)])
Ess2[4] = len(df2[(df2['TM-score'] >= 0.8) & (df2['TM-score'] < 0.9) & (df2['Num_ess'] > 0)])
Ess2[5] = len(df2[(df2['TM-score'] >= 0.9) & (df2['Num_ess'] > 0)])


df_2 = pd.DataFrame({'TM-score': ['0.474-0.5', '0.5-0.6', '0.6-0.7', '0.7-0.8', '0.8-0.9', '0.9-1'], 
                    '< 10 intersections': OK.tolist(), 
                    '>= 10 intersectios ': No_Ess.tolist(), 
                    'Has essential intersections': Ess.tolist()})

df_22 = pd.DataFrame({'TM-score': ['0.474-0.5', '0.5-0.6', '0.6-0.7', '0.7-0.8', '0.8-0.9', '0.9-1'], 
                    '< 10 intersections': OK2.tolist(), 
                    '>= 10 intersectios ': No_Ess2.tolist(), 
                    'Has essential intersections': Ess2.tolist()})

# Add a new column to each DataFrame to indicate the DataFrame it belongs to
df_2['DataFrame'] = 'Max length: 18'
df_22['DataFrame'] = 'Max length: 50'

# Concatenate the two DataFrames
df_concat = pd.concat([df_2, df_22])

# Reshape the DataFrame
df_concat_melt = df_concat.melt(id_vars=['TM-score', 'DataFrame'], var_name='Error measure', value_name='Percentage of models')

# Create the bar chart
fig = px.bar(df_concat_melt, x='TM-score', y='Percentage of models', color='Error measure', barmode='stack', 
             facet_col='DataFrame', color_discrete_sequence=["#77dd77", "#6699cc", "#ff6961"])

# Update layout to make the plot wider
fig.update_layout(autosize=False, width=1200)

fig.show()


# --------- Now percentage --------------------
# Reshape the DataFrame
df_concat_melt = df_concat.melt(id_vars=['TM-score', 'DataFrame'], var_name='Error measure', value_name='Percentage of models')

# Calculate the total for each 'TM-score' category
total = df_concat_melt.groupby(['TM-score', 'DataFrame'])['Percentage of models'].transform('sum')

# Convert 'Number of models' to percentage
df_concat_melt['Percentage of models'] = df_concat_melt['Percentage of models'] / total * 100

# Create the bar chart
fig = px.bar(df_concat_melt, x='TM-score', y='Percentage of models', color='Error measure', facet_col='DataFrame', barmode='stack', 
             color_discrete_sequence=["#77dd77", "#6699cc", "#ff6961"])

# Update layout to make the plot wider
fig.update_layout(autosize=False, width=1200)

fig.show()




#fig = px.bar(df, x='TM-score', y=df.columns[1:], barmode='stack')

# Reshape the DataFrame
#df_2_melt = df_2.melt(id_vars='TM-score', var_name='Error measure', value_name='Percentage of models')

# Create the bar chart
#fig = px.bar(df_2_melt, x='TM-score', y='Percentage of models', color='Error measure', barmode='stack', 
#             color_discrete_sequence=["#77dd77", "#6699cc", "#ff6961"])

# Update layout to make the plot wider
#fig.update_layout(autosize=False, width=1200)

#fig.show()

# --------- Now procentage --------------------
# Reshape the DataFrame
#df2_melt = df_2.melt(id_vars='TM-score', var_name='Error measure', value_name='Percentage of models')

# Calculate the total for each 'TM-score' category
#total = df2_melt.groupby('TM-score')['Percentage of models'].transform('sum')

# Convert 'Number of models' to percentage
#df2_melt['Percentage of models'] = df2_melt['Percentage of models'] / total * 100

# Create the bar chart
#fig = px.bar(df2_melt, x='TM-score', y='Percentage of models', color='Error measure', barmode='stack', 
#             color_discrete_sequence=["#77dd77", "#6699cc", "#ff6961"])

# Update layout to make the plot wider
#fig.update_layout(autosize=False, width=1200)

#fig.show()
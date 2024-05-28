from CifFile import ReadCif

# Read the CIF file
cif_file = ReadCif('C:/Users/Kapta/Documents/Skole/DTU/6.semester/BP/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/Alphafold/fold_2024_05_28_10_03_model_0.cif')

# Extract the data block (replace 'sample' with your actual data block name if different)
data_block = cif_file.first_block()

# Extract PDB-related atom site information
group_pdb = data_block['_atom_site.group_PDB']
atom_id = data_block['_atom_site.id']
type_symbol = data_block['_atom_site.type_symbol']
label_atom_id = data_block['_atom_site.label_atom_id']
label_alt_id = data_block['_atom_site.label_alt_id']
label_comp_id = data_block['_atom_site.label_comp_id']
label_asym_id = data_block['_atom_site.label_asym_id']
label_entity_id = data_block['_atom_site.label_entity_id']
label_seq_id = data_block['_atom_site.label_seq_id']
pdbx_PDB_ins_code = data_block['_atom_site.pdbx_PDB_ins_code']
cartn_x = data_block['_atom_site.Cartn_x']
cartn_y = data_block['_atom_site.Cartn_y']
cartn_z = data_block['_atom_site.Cartn_z']
occupancy = data_block['_atom_site.occupancy']
b_iso_or_equiv = data_block['_atom_site.B_iso_or_equiv']
auth_seq_id = data_block['_atom_site.auth_seq_id']
auth_asym_id = data_block['_atom_site.auth_asym_id']
pdbx_PDB_model_num = data_block['_atom_site.pdbx_PDB_model_num']

# Write to PDB format
pdb_lines = []
for i in range(len(group_pdb)):
    pdb_line = (
        f"{group_pdb[i]:<6}" +
        f"{int(atom_id[i]):>5} " +
        f"{label_atom_id[i]:<4}{label_alt_id[i]:1}" +
        f"{label_comp_id[i]:>3} " +
        f"{label_asym_id[i]:1}" +
        f"{int(label_entity_id[i]):>4}" +
        f"{pdbx_PDB_ins_code[i]:1}   " +
        f"{float(cartn_x[i]):>8.3f}" +
        f"{float(cartn_y[i]):>8.3f}" +
        f"{float(cartn_z[i]):>8.3f}" +
        f"{float(occupancy[i]):>6.2f}" +
        f"{float(b_iso_or_equiv[i]):>6.2f}  " +
        f"{int(auth_seq_id[i]):>4}  " +
        f"{label_asym_id[i]:1} " +
        f"{int(pdbx_PDB_model_num[i]):>2}\n"
    )
    pdb_lines.append(pdb_line)

# Save to PDB file
pdb_content = ''.join(pdb_lines)
with open('output.pdb', 'w') as pdb_file:
    pdb_file.write(pdb_content)

print("PDB file saved as 'output.pdb'")

from Bio.PDB import PDBParser
import pandas as pd
from Bio.PDB.Polypeptide import PPBuilder

def one_PDB_to_seq(PDB_filename):

    p = PDBParser(QUIET=True)
    s1 = p.get_structure("s1", PDB_filename)

    ca1 = [(atom.full_id[2], atom.full_id[3][1], *atom.get_coord()) for atom in s1.get_atoms() if atom.full_id[4][0] == "CA"]
    df1 = pd.DataFrame(ca1, columns=["chain", "residue_number", "x", "y", "z"])
    

    df1["P"] = df1[["x", "y", "z"]].values.tolist()

    P1 = {}
    
    for chain in df1["chain"].unique():
        P1["Chain_" + chain] = df1[df1["chain"] == chain]["P"]

    # Get the sequence
    seq1 = {}
    ppb=PPBuilder()
    i = 0
    for pp in ppb.build_peptides(s1):
        seq1["Chain_" + df1["chain"].unique()[i]] = pp.get_sequence()
        i += 1
    i = 0
    return P1, seq1, s1


def two_PDB_to_seq(PDB1_filename, PDB2_filename):

    p = PDBParser(QUIET=True)
    s1 = p.get_structure("s1", PDB1_filename)
    s2 = p.get_structure("s2", PDB2_filename)

    ca1 = [(atom.full_id[2], atom.full_id[3][1], *atom.get_coord()) for atom in s1.get_atoms() if atom.full_id[4][0] == "CA"]
    df1 = pd.DataFrame(ca1, columns=["chain", "residue_number", "x", "y", "z"])
    
    ca2 = [(atom.full_id[2], atom.full_id[3][1], *atom.get_coord()) for atom in s2.get_atoms() if atom.full_id[4][0] == "CA"]
    df2 = pd.DataFrame(ca2, columns=["chain", "residue_number", "x", "y", "z"])

    x_mean = df1["x"].mean()
    y_mean = df1["y"].mean()
    z_mean = df1["z"].mean()

    df1["x"] = df1["x"] - x_mean
    df1["y"] = df1["y"] - y_mean
    df1["z"] = df1["z"] - z_mean

    x_mean = df2["x"].mean()
    y_mean = df2["y"].mean()
    z_mean = df2["z"].mean()

    df2["x"] = df2["x"] - x_mean
    df2["y"] = df2["y"] - y_mean
    df2["z"] = df2["z"] - z_mean

    df1["P"] = df1[["x", "y", "z"]].values.tolist()
    df2["P"] = df2[["x", "y", "z"]].values.tolist()

    P1 = {}
    P2 = {}
    if len(df1["chain"].unique()) != len(df2["chain"].unique()):
        print("The number of chains is different between the two structures")
        return
    
    for chain in df1["chain"].unique():
        P1["Chain_" + chain] = df1[df1["chain"] == chain]["P"]
        P2["Chain_" + chain] = df2[df2["chain"] == chain]["P"]

    # Get the sequence
    seq1 = {}
    seq2 = {}
    ppb=PPBuilder()
    i = 0
    for pp in ppb.build_peptides(s1):
        seq1["Chain_" + df1["chain"].unique()[i]] = pp.get_sequence()
        i += 1
    i = 0
    for pp in ppb.build_peptides(s2):
        seq2["Chain_" + df2["chain"].unique()[i]] = pp.get_sequence()
        i += 1

    return P1, P2, seq1, seq2, s1, s2
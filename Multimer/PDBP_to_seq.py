from Bio.PDB import PDBParser
import pandas as pd
from Bio.Seq import Seq
from Bio.PDB.Polypeptide import PPBuilder, CaPPBuilder
import numpy as np

def one_PDB_to_seq(PDB_filename):

    p = PDBParser(QUIET=True)
    s1 = p.get_structure("s1", PDB_filename)

    ca1 = [(atom.full_id[2], atom.full_id[3][1], *atom.get_coord()) for atom in s1.get_atoms() if atom.full_id[4][0] == "CA"]
    df1 = pd.DataFrame(ca1, columns=["chain", "residue_number", "x", "y", "z"])
    

    # df1["x"] = df1["x"] - s1.center_of_mass()[0]
    # df1["y"] = df1["y"] - s1.center_of_mass()[1]
    # df1["z"] = df1["z"] - s1.center_of_mass()[2]

    
    df2 = df1
    df1["P"] = df1[["x", "y", "z"]].values.tolist()

    P1 = {}
    chain_com = {}
    for chain in df1["chain"].unique():
        P1["Chain_" + chain] = df1[df1["chain"] == chain]["P"]
        chain_data_x = df1[df1["chain"] == chain]["x"].mean()
        chain_data_y = df1[df1["chain"] == chain]["y"].mean()
        chain_data_z = df1[df1["chain"] == chain]["y"].mean()
        chain_com["Chain_" + chain] = [chain_data_x, chain_data_y, chain_data_z]

    # Get the sequence
    seq1 = {}
    tot_seq1 = Seq("")
    ppb=PPBuilder()
    if len(ppb.build_peptides(s1)) != len(df1["chain"].unique()):
        ppb=CaPPBuilder()
        

    i = 0
    for pp in ppb.build_peptides(s1):
        seq1["Chain_" + df1["chain"].unique()[i]] = pp.get_sequence()
        tot_seq1 += pp.get_sequence()
        #CaPP
        i += 1
    
    return P1, seq1, s1, tot_seq1, chain_com


def two_PDB_to_seq(PDB1_filename, PDB2_filename):

    P1, seq1, s1, tot_seq1, chain_com1 = one_PDB_to_seq(PDB1_filename)
    P2, seq2, s2, tot_seq2, chain_com2 = one_PDB_to_seq(PDB2_filename)


    if len(seq1) != len(seq2):
        print("The number of chains is different between the two structures")

    return P1, P2, seq1, seq2, s1, s2, tot_seq1, tot_seq2, chain_com1, chain_com2
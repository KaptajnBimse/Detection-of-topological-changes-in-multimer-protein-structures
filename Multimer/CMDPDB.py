from Bio import PDB
from Bio.SeqUtils import seq1

class ChainSequence:
    def __init__(self, chain_id, sequence):
        self.chain_id = chain_id
        self.sequence = sequence

    def __repr__(self):
        return f"Chain {self.chain_id}: {self.sequence}"

# Function to extract sequences from all chains and store in ChainSequence objects
def extract_sequences_from_pdb(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    chain_sequences = []
    
    for model in structure:
        for chain in model:
            seq = ""
            for residue in chain:
                if PDB.is_aa(residue, standard=True):
                    seq += seq1(residue.resname)
            chain_sequences.append(ChainSequence(chain.id, seq))
    
    return chain_sequences

# Path to your PDB file
pdb_file_path = "/Users/agb/Desktop/Bachelorprojekt/Detection-of-topological-changes-in-multimer-protein-structures/Multimer/examples/Multimer PDB/fixed.pdb"  # Replace with your PDB file path

# Extract sequences and store them in ChainSequence objects
chain_sequences = extract_sequences_from_pdb(pdb_file_path)

# Print the sequences for each chain
for chain_seq in chain_sequences:
    print(chain_seq)


import os
from Bio import PDB
from Bio.PDB.Polypeptide import is_aa
from collections import defaultdict
import numpy as np

def classify_residue(residue):
    polar_positive_residues = ["ARG", "LYS", "HIS"]
    polar_negative_residues = ["GLU", "ASP"]
    
    if residue.get_resname() in polar_positive_residues:
        return "polar positive"
    elif residue.get_resname() in polar_negative_residues:
        return "polar negative"
    else:
        return "neutral"

def process_protein(protein, n1_cumulative_counts, n5_cumulative_counts):
    fad_candidates = []  

    for model in protein:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == "FMN":
                    n1_coords = None
                    n5_coords = None

                    for atom in residue:
                        if atom.get_name() == "N1":
                            n1_coords = atom.get_coord()
                        elif atom.get_name() == "N5":
                            n5_coords = atom.get_coord()

                    if n1_coords is not None and n5_coords is not None:
                        print(f"N1 Coordinates for bifurcator {residue.get_full_id()}:", n1_coords)
                        print(f"N5 Coordinates for bifurcator {residue.get_full_id()}:", n5_coords)

                        counted_residues_n1 = set()  # Keep track of counted residues for N1
                        counted_residues_n5 = set()  # Keep track of counted residues for N5

                        for distance in range(4, 13):  # Loop through distances from 4 to 12 Angstroms
                            for atom in protein.get_atoms():
                                if is_aa(atom.get_parent()):
                                    dist_n1 = np.linalg.norm(atom.get_coord() - n1_coords)
                                    dist_n5 = np.linalg.norm(atom.get_coord() - n5_coords)

                                    if dist_n1 <= distance and atom.get_parent().id not in counted_residues_n1:
                                        classification = classify_residue(atom.get_parent())
                                        if classification == "polar positive":
                                            n1_cumulative_counts[distance][0] += 1
                                        elif classification == "polar negative":
                                            n1_cumulative_counts[distance][1] += 1
                                        else:
                                            n1_cumulative_counts[distance][2] += 1
                                        counted_residues_n1.add(atom.get_parent().id)

                                    if dist_n5 <= distance and atom.get_parent().id not in counted_residues_n5:
                                        classification = classify_residue(atom.get_parent())
                                        if classification == "polar positive":
                                            n5_cumulative_counts[distance][0] += 1
                                        elif classification == "polar negative":
                                            n5_cumulative_counts[distance][1] += 1
                                        else:
                                            n5_cumulative_counts[distance][2] += 1
                                        counted_residues_n5.add(atom.get_parent().id)

# Assuming you have a directory called 'proteins_to_use' containing PDB files
protein_files = [f for f in os.listdir('proteins_to_use') if f.endswith('.pdb')]

n1_cumulative_counts = defaultdict(lambda: [0, 0, 0])  # [Polar Positive, Polar Negative, Neutral]
n5_cumulative_counts = defaultdict(lambda: [0, 0, 0])  # [Polar Positive, Polar Negative, Neutral]

# Load the protein structures
for protein_file in protein_files:
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(protein_file, os.path.join('proteins_to_use', protein_file))
    process_protein(structure, n1_cumulative_counts, n5_cumulative_counts)

# Write the results to files
def write_counts_to_file(filename, cumulative_counts):
    with open(filename, 'w') as counts_file:
        counts_file.write("Distance,Polar Positive,Polar Negative,Neutral\n")
        for distance, counts in cumulative_counts.items():
            counts_file.write(f"{distance},{counts[0]},{counts[1]},{counts[2]}\n")

write_counts_to_file('n1_counts.txt', n1_cumulative_counts)
write_counts_to_file('n5_counts.txt', n5_cumulative_counts)


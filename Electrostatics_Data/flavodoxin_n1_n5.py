import os
from Bio import PDB
from Bio.PDB.Polypeptide import is_aa
import numpy as np

def classify_residue(residue):
    nonpolar_residues = ["ALA", "VAL", "LEU", "ILE", "MET", "PRO", "PHE", "TRP", "TYR", "GLY", "CYS", "MSE"]
    
    if residue.get_resname() in nonpolar_residues:
        return "nonpolar"
    else:
        return "polar"

def process_protein(protein):
    fmn_n1_coords = None
    fmn_n5_coords = None

    for atom in protein.get_atoms():
        if atom.get_parent().get_resname() == "FMN":
            if atom.get_name() == "N1":
                fmn_n1_coords = atom.get_coord()
            elif atom.get_name() == "N5":
                fmn_n5_coords = atom.get_coord()

    if fmn_n1_coords is not None and fmn_n5_coords is not None:
        distances_n1 = {}
        distances_n5 = {}

        for distance in range(4, 13):
            for atom in protein.get_atoms():
                if is_aa(atom.get_parent()):
                    dist_n1 = np.linalg.norm(atom.get_coord() - fmn_n1_coords)
                    if dist_n1 <= distance:
                        if distance not in distances_n1:
                            distances_n1[distance] = []
                        distances_n1[distance].append((atom.get_parent().get_resname(), dist_n1))

                    dist_n5 = np.linalg.norm(atom.get_coord() - fmn_n5_coords)
                    if dist_n5 <= distance:
                        if distance not in distances_n5:
                            distances_n5[distance] = []
                        distances_n5[distance].append((atom.get_parent().get_resname(), dist_n5))

        return distances_n1, distances_n5
    else:
        return None, None

# Assuming you have a directory called 'proteins_to_use' containing PDB files
protein_files = [f for f in os.listdir('proteins_to_use') if f.endswith('.pdb')]

# Load the protein structures
flavodoxins = []
for protein_file in protein_files:
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(protein_file, os.path.join('proteins_to_use', protein_file))
    flavodoxins.append(structure)

# Process each protein structure
for protein in flavodoxins:
    n1_distances, n5_distances = process_protein(protein)
    if n1_distances is not None and n5_distances is not None:
        for distance in range(4, 13):
            if distance in n1_distances:
                with open(f'n1_distances_{distance}A.txt', 'a') as n1_file:
                    for resname, dist in n1_distances[distance]:
                        n1_file.write(f"{protein.id} {resname} {dist:.2f}\n")

            if distance in n5_distances:
                with open(f'n5_distances_{distance}A.txt', 'a') as n5_file:
                    for resname, dist in n5_distances[distance]:
                        n5_file.write(f"{protein.id} {resname} {dist:.2f}\n")


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

def process_bifurcator(protein):
    fad_coords = None

    # Identify the bifurcating FAD with an HN5 atom
    for atom in protein.get_atoms():
        if atom.get_parent().get_resname() == "FAD":
            if atom.get_name() == "HN5":
                fad_coords = atom.get_coord()
                break

    if fad_coords is not None:
        distances_n1 = {}
        distances_n5 = {}

        for distance in range(4, 13):
            for atom in protein.get_atoms():
                if is_aa(atom.get_parent()):
                    dist_n1 = np.linalg.norm(atom.get_coord() - fad_coords)
                    if dist_n1 <= distance:
                        if distance not in distances_n1:
                            distances_n1[distance] = []
                        distances_n1[distance].append((atom.get_parent().get_resname(), dist_n1))

                    dist_n5 = np.linalg.norm(atom.get_coord() - fad_coords)
                    if dist_n5 <= distance:
                        if distance not in distances_n5:
                            distances_n5[distance] = []
                        distances_n5[distance].append((atom.get_parent().get_resname(), dist_n5))

        return distances_n1, distances_n5
    else:
        return None, None

# Assuming you have a directory called 'bifurcators_to_use' containing PDB files
bifurcator_files = [f for f in os.listdir('bifurcators_to_use') if f.endswith('.pdb')]

# Load the bifurcator structures
bifurcators = []
for bifurcator_file in bifurcator_files:
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(bifurcator_file, os.path.join('bifurcators_to_use', bifurcator_file))
    bifurcators.append(structure)

# Process each bifurcator structure
for bifurcator in bifurcators:
    n1_distances, n5_distances = process_bifurcator(bifurcator)
    if n1_distances is not None and n5_distances is not None:
        for distance in range(4, 13):
            if distance in n1_distances:
                with open(f'n1_distances_{distance}A.txt', 'a') as n1_file:
                    for resname, dist in n1_distances[distance]:
                        n1_file.write(f"{bifurcator.id} {resname} {dist:.2f}\n")

            if distance in n5_distances:
                with open(f'n5_distances_{distance}A.txt', 'a') as n5_file:
                    for resname, dist in n5_distances[distance]:
                        n5_file.write(f"{bifurcator.id} {resname} {dist:.2f}\n")


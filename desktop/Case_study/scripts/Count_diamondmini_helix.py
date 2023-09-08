import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

# Subdirectories
Subdirectory = ["e6", "e7", "e8", "e9", "e10", "e11"]



# Format the ratio as "50:50" style
def format_ratio(rratio, lratio):
    return f"{int(rratio * 100)}:{int(lratio * 100)}"

# Calculate dihedral angles from an sdf file
def calculate_dihedrals(sdf_file):
    suppl = Chem.SDMolSupplier(sdf_file)
    dihedral_angles_all = []

    for mol in suppl:
        if mol is not None:
            conf = mol.GetConformer(0)
            num_atoms = mol.GetNumAtoms()
            for i in range(num_atoms - 3):
                for j in range(i + 1, num_atoms - 2):
                    for k in range(j + 1, num_atoms - 1):
                        for l in range(k + 1, num_atoms):
                            dihedral = Chem.rdMolTransforms.GetDihedralDeg(conf, i, j, k, l)
                            dihedral_angles_all.append(dihedral)
    return dihedral_angles_all

# Calculate helix ratios from dihedral angles with a tolerance range
def find_helix_ratio(dihedral_angles, sequence_length=4, tolerance=5):
    left_helix_pattern = [180, 60, 180, 60]
    right_helix_pattern = [-60, 180, -60, 180]

    def is_helix_match(angle_sequence, pattern, tolerance):
        return all(abs(a - p) < tolerance for a, p in zip(angle_sequence, pattern))

    # Creating subsequences of the length of the pattern
    subsequences = [dihedral_angles[i:i + sequence_length] for i in range(len(dihedral_angles) - sequence_length + 1)]

    left_helix_count = sum(is_helix_match(subseq, left_helix_pattern, tolerance) for subseq in subsequences)
    right_helix_count = sum(is_helix_match(subseq, right_helix_pattern, tolerance) for subseq in subsequences)

    total_count = len(subsequences)

    left_ratio = left_helix_count / total_count if total_count > 0 else 0
    right_ratio = right_helix_count / total_count if total_count > 0 else 0

    return left_ratio, right_ratio


# Initialize a list to store molecule data
molecule_data = []

# Iterate through subdirectories
for subdir in Subdirectory:
    dihedral_angles_all = []
    for filename in os.listdir(subdir):
        if filename.endswith(".sdf"):
            sdf_file = os.path.join(subdir, filename)
            dihedral_angles = calculate_dihedrals(sdf_file)
            dihedral_angles_all.extend(dihedral_angles)

    left_ratio, right_ratio = find_helix_ratio(dihedral_angles_all)

    total_helix_ratio = left_ratio + right_ratio
    formatted_ratio = "Undefined"  # Default value

    if total_helix_ratio > 0:
        formatted_ratio = format_ratio(right_ratio / total_helix_ratio, left_ratio / total_helix_ratio)

    molecule_data.append({
        "Molecular Name": subdir,
        "Total Helix Ratio": total_helix_ratio,
        "Left Helix Ratio": left_ratio,
        "Right Helix Ratio": right_ratio,
        "Formatted Helix Ratio": formatted_ratio
    })

# Create DataFrame and save to CSV
columns_order = ["Molecular Name", "Total Helix Ratio", "Left Helix Ratio", "Right Helix Ratio", "Formatted Helix Ratio"]
molecule_df = pd.DataFrame(molecule_data, columns=columns_order)
molecule_df.to_csv("molecule_helix_ratios.csv", index=False)


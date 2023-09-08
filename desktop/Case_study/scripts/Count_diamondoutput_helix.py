import os
import csv
import math

# Function to find helix
def find_helix(data):
    left_helix_pattern = [180, 60, 180, 60]  # Left-handed helix (M): t g+ t g+ t
    right_helix_pattern = [-60, 180, -60, 180]  # Right-handed helix: g− t g− t

    left_helix_count = 0
    right_helix_count = 0
    total_count = len(data)

    for dihedral in data:
        for i in range(len(dihedral) - len(left_helix_pattern) + 1):
            if dihedral[i:i + len(left_helix_pattern)] == left_helix_pattern:
                left_helix_count += 1
            if dihedral[i:i + len(right_helix_pattern)] == right_helix_pattern:
                right_helix_count += 1

    left_helix_ratio = left_helix_count / total_count if total_count > 0 else 0
    right_helix_ratio = right_helix_count / total_count if total_count > 0 else 0

    return left_helix_ratio, right_helix_ratio, left_helix_count, right_helix_count


# Function to format the ratio as "50:50" style
def format_ratio(rratio, lratio):
    return f"{int(rratio * 100)}:{int(lratio * 100)}"


# Subdirectories
Subdirectory = ["e6", "e7", "e8", "e9", "e10", "e11"]

# Prepare a CSV file to store results
with open("helix_ratios.csv", "w", newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(
        ["Molecular Name", "Total Helix Ratio", "Left Helix Ratio", "Right Helix Ratio", "Formatted Helix Ratio"])

    for subdir in Subdirectory:
        filepath = os.path.join(subdir, "running_record.log")

        # Initialize data
        data = []

        # Read each file
        try:
            with open(filepath, "r") as f:
                lines = f.readlines()
                line_num = 0
                for line in lines:
                    line_num += 1
                    if "Conformation, Energy" in line:
                        parts = line.split(']')
                        if len(parts) < 3:
                            print(f"Insufficient parts in line {line_num}.")
                            continue
                        dihedral_info = [float(x) for x in parts[1].strip()[1:].split(",")]
                        print("dihedral_info:", dihedral_info)
                        data.append(dihedral_info)

            # Find helix ratio
            left_ratio, right_ratio, left_count, right_count = find_helix(data)
            total_ratio = left_ratio + right_ratio
            # Write to CSV
            csvwriter.writerow([subdir, total_ratio, left_ratio, right_ratio, format_ratio(right_ratio / total_ratio, left_ratio / total_ratio)])

            print(f"Subdirectory {subdir}:")
            print(f"Total Helix Ratio: {total_ratio}")
            print(f"Left Helix Ratio: {left_ratio}")
            print(f"Right Helix Ratio: {right_ratio}")
            print(f"Formatted Helix Ratio: {format_ratio(right_ratio / total_ratio, left_ratio / total_ratio)}")
            print(f"Left Helix Count: {left_count}")
            print(f"Right Helix Count: {right_count}\n")

        except FileNotFoundError:
            print(f"File not found for subdirectory {subdir}")

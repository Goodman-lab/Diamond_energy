import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import Draw
from rdkit.Chem.rdmolfiles import SDWriter
from rdkit.Chem import TorsionFingerprints
from rdkit.Chem.Draw import rdDepictor
from rdkit import rdBase
import os
import pandas as pd
from datetime import datetime
import re
import matplotlib.pyplot as plt


# The molecular name in acyclic alkane test set
# One can change the following list to target test set
name=["c4", "c5", "c5n", "c6", "c6n", "c7", "c7c", "c7n", "c8", "c8a",
      "c8s", "c8si", "c9", "c9n", "c9s", "c10n", "c12p", "c12p0", "c12p1", "c15",
      "c15helix", "c15a", "c15b", "c15c", "c15d", "c15e", "c15f", "c15g", "c16", "c17",
      "c17a", "c17l", "c18", "c18n", "c19", "c20n", "c21", "c26"]
total_number=len(name)
# The comparing experiment names of Macromodel with algorithm and force filed
Maestro_test=["MCMM_OPLS4", "MCMM_MM3", "SPMC_OPLS4", "SPMC_MM3",
              "LMCS_OPLS4", "LMCS_MM3", "MCLMCS_OPLS4", "MCLMCS_MM3"]
total_method=len(Maestro_test)

# Function to calculate running time from log lines
def read_eachMaestrotest():
    # Add the calculate_time function
    def calculate_time(start_line, end_line):
        start_time_str = re.search(r"Starting Time (.+)$", start_line).group(1)
        end_time_str = re.search(r"termination\s+(.+)$", end_line).group(1)

        start_time = datetime.strptime(start_time_str, '%d-%b-%Y %H:%M:%S')
        end_time = datetime.strptime(end_time_str, '%d-%b-%Y %H:%M:%S')

        return (end_time - start_time).seconds

    for test_method in Maestro_test:
        final_df = pd.DataFrame(columns=['Molecule', 'Lowest_Energy', 'Unique_Structures', 'Total_Running_Time', 'Average_Running_Time'])

        for molecule in name:
            file_path = os.path.join(".", test_method, f"{molecule}_{test_method}.log")
            print("file_path:", file_path)
            try:
                with open(file_path, "r", encoding='utf-8') as f:
                    lines = f.readlines()
                    Lowest_energy_mae = None
                    Output_number = None
                    start_line = None
                    end_line = None

                    for line in lines:
                        line = line.strip()
                        # Debug: Print the line being processed
                        #print("Processing line: ", line.strip())
                        # Using the regular expression to search for the pattern in each line
                        conformation_match = re.search(
                            r'Conformation\s+1\s+\(\s+([\d.]+)\s+kJ/mol\)\s+was\s+found\s+[\d]+\s+times', line)

                        if conformation_match:
                            # Capturing the energy value using .group(1)
                            try:
                                Lowest_energy_mae = float(conformation_match.group(1))
                                print(f"Lowest_energy_mae successfully captured: {Lowest_energy_mae}")
                            except ValueError as e:
                                print(f"Skipping line due to conversion error: {e}")
                        elif "minimized with good convergence" in line:
                            try:
                                Output_number=int(line.split("minimized")[0].strip())
                                print("Output_number:", Output_number)
                            except ValueError as e:
                                print(f"Skipping line due to conversion error: {e}")
                                continue
                        elif "BatchMin: normal termination" in line:
                            end_line = line
                            print("end_line:", end_line)
                        elif "Starting Time" in line:
                            start_line = line
                            print("start_line:", start_line)

                    running_time = calculate_time(start_line, end_line) if start_line and end_line else None
                    average_running_time = running_time / Output_number if Output_number and running_time else None

                    new_row = pd.DataFrame({
                        'Molecule': [molecule],
                        'Lowest_Energy': [Lowest_energy_mae],
                        'Unique_Structures': [Output_number],
                        'Total_Running_Time': [running_time],
                        'Average_Running_Time': [average_running_time]
                    })

                    final_df = pd.concat([final_df, new_row], ignore_index=True)
                f.close()
            except FileNotFoundError:
                print(f"File {file_path} not found.")
                continue
        print("final_df:", final_df)
        # Save the DataFrame as a CSV file in the corresponding subdirectory
        final_df.to_csv(os.path.join(".", test_method, 'results.csv'), index=False)
        break


# Run the function
#read_eachMaestrotest()


#split_directory()
#run_maestroscript_parallel()

def draw_results_initial():
    import pandas as pd
    import matplotlib.pyplot as plt
    import os
    def add_jitter(arr, factor=0.005):
        """
        Adds a small jitter to the values in arr.
        """
        return arr + np.random.randn(len(arr)) * factor

    # List of filenames and corresponding plot properties
    files_and_properties = [
        {"filename": "Diamond_OPLS4.csv", "color": "red", "label": "Diamond_OPLS4", "marker": 'o', "markersize": 3},
        {"filename": "MCMM_OPLS4.csv", "color": "blue", "label": "MCMM_OPLS4", "marker": 's', "markersize": 9},
        {"filename": "MCLMCS_OPLS4.csv", "color": "yellow", "label": "MCLMCS_OPLS4", "marker": '^', "markersize": 7},
        {"filename": "SPMC_OPLS4.csv", "color": "purple", "label": "SPMC_OPLS4", "marker": 'D', "markersize": 5},
        {"filename": "LMCS_OPLS4.csv", "color": "green", "label": "LMCS_OPLS4", "marker": 'x', "markersize": 3}
    ]

    # Initialize an empty dataframe to hold the common molecule names across all files
    common_molecule_names = None

    # Read each file and plot
    for prop in files_and_properties:
        filename = prop["filename"]
        color = prop["color"]
        label = prop["label"]

        if os.path.exists(filename):
            df = pd.read_csv(filename)
            # Convert the second column to numeric, invalid parsing will be set as NaN
            df.iloc[:, 1] = pd.to_numeric(df.iloc[:, 1], errors='coerce')

            # Drop rows where any of the two required columns has NaN
            df.dropna(subset=[df.columns[0], df.columns[1]], inplace=True)

            # If common_molecule_names is still None, initialize it
            if common_molecule_names is None:
                common_molecule_names = set(df.iloc[:, 0])
            else:
                # Update common_molecule_names to only contain names common to all files read so far
                common_molecule_names = common_molecule_names.intersection(set(df.iloc[:, 0]))

    # Filter and plot data only for common molecule names
    for prop in files_and_properties:
        filename = prop["filename"]
        color = prop["color"]
        label = prop["label"]

        if os.path.exists(filename):
            df = pd.read_csv(filename)
            df.iloc[:, 1] = pd.to_numeric(df.iloc[:, 1], errors='coerce')
            df = df[df.iloc[:, 0].isin(common_molecule_names)]

            plt.plot(df.iloc[:, 0], df.iloc[:, 1], linestyle='-', marker='o', color=color, label=label)

    # Configure plot
    plt.title('Diamond vs other methods in lowest conformer found')
    plt.xlabel('Molecule Name')
    plt.ylabel('Lowest Energy')
    plt.legend()

    # Show the plot
    plt.show()


import matplotlib.pyplot as plt
import numpy as np  # For adding jitter

def draw_results():
    def add_jitter(arr, factor=0.005):
        """
        Adds a small jitter to the values in arr.
        """
        return arr + np.random.randn(len(arr)) * factor

    # List of filenames and corresponding plot properties
    files_and_properties = [
        {"filename": "Diamond_OPLS4.csv", "color": "red", "label": "Diamond_OPLS4", "marker": 'o', "markersize": 6},
        {"filename": "MCMM_OPLS4.csv", "color": "blue", "label": "MCMM_OPLS4", "marker": 's', "markersize": 8},
        {"filename": "MCLMCS_OPLS4.csv", "color": "yellow", "label": "MCLMCS_OPLS4", "marker": '^',"markersize": 8},
        {"filename": "SPMC_OPLS4.csv", "color": "purple", "label": "SPMC_OPLS4", "marker": 'D', "markersize": 8},
        {"filename": "LMCS_OPLS4.csv", "color": "green", "label": "LMCS_OPLS4", "marker": 'x', "markersize": 8}
    ]

    # Initialize an empty dataframe to hold the common molecule names across all files
    common_molecule_names = None

    # Read each file and plot
    for prop in files_and_properties:
        filename = prop["filename"]

        if os.path.exists(filename):
            df = pd.read_csv(filename)
            df.iloc[:, 1] = pd.to_numeric(df.iloc[:, 1], errors='coerce')
            df.dropna(subset=[df.columns[0], df.columns[1]], inplace=True)

            if common_molecule_names is None:
                common_molecule_names = set(df.iloc[:, 0])
            else:
                common_molecule_names = common_molecule_names.intersection(set(df.iloc[:, 0]))

    # Filter and plot data only for common molecule names
    for prop in files_and_properties:
        filename = prop["filename"]
        color = prop["color"]
        label = prop["label"]
        marker = prop["marker"]
        markersize = prop["markersize"]

        if os.path.exists(filename):
            df = pd.read_csv(filename)
            df.iloc[:, 1] = pd.to_numeric(df.iloc[:, 1], errors='coerce')
            df = df[df.iloc[:, 0].isin(common_molecule_names)]

            # Sort the dataframe by the 'Molecule' column to have a consistent x-axis
            df.sort_values(by=df.columns[0], inplace=True)

            y_values = df.iloc[:, 1]

            # Add a small jitter to the y-values for all lines except "Diamond_OPLS4"
            if label != "Diamond_OPLS4":
                y_values = add_jitter(y_values)

            plt.plot(df.iloc[:, 0], y_values, linestyle='-', marker=marker, color=color, markersize=markersize,
                     label=label)

    # Configure plot
    plt.title('Diamond vs other methods in lowest conformer found')
    plt.xlabel('Molecule Name')
    plt.ylabel('Lowest Energy')
    plt.legend()

    # Show the plot
    plt.show()



draw_results()


def draw_results_test():
    import pandas as pd
    import matplotlib.pyplot as plt
    import os

    # List of filenames and corresponding plot properties
    files_and_properties = [
        {"filename": "Diamond_OPLS4.csv", "color": "red", "label": "Diamond_OPLS4", "marker": 'o', "markersize": 6},
        {"filename": "MCMM_OPLS4.csv", "color": "blue", "label": "MCMM_OPLS4", "marker": 's', "markersize": 8},
        {"filename": "MCLMCS_OPLS4.csv", "color": "yellow", "label": "MCLMCS_OPLS4", "marker": '^', "markersize": 8},
        {"filename": "SPMC_OPLS4.csv", "color": "purple", "label": "SPMC_OPLS4", "marker": 'D', "markersize": 8},
        {"filename": "LMCS_OPLS4.csv", "color": "green", "label": "LMCS_OPLS4", "marker": 'x', "markersize": 8}
    ]

    # Initialize an empty dataframe to hold the common molecule names across all files
    common_molecule_names = None
    if common_molecule_names is None:
        print("No common molecule names found. Please check if the files are correct.")
        return

    # Filter and plot data only for common molecule names
    for prop in files_and_properties:
        filename = prop["filename"]
        color = prop["color"]
        label = prop["label"]

        if os.path.exists(filename):
            df = pd.read_csv(filename)
            df.iloc[:, 1] = pd.to_numeric(df.iloc[:, 1], errors='coerce')

            # Only keep rows with 'Molecule' names that are in the common set
            df = df[df.iloc[:, 0].isin(common_molecule_names)]

            # Sort the dataframe by the 'Molecule' column
            df.sort_values(by=df.columns[0], inplace=True)

            plt.plot(df.iloc[:, 0], df.iloc[:, 1], linestyle='-', marker='o', color=color, label=label)

    # Configure plot
    plt.title('Diamond vs other methods in lowest conformer found')
    plt.xlabel('Molecule Name')
    plt.ylabel('Lowest Energy')
    plt.legend()

    # Show the plot
    plt.show()


#draw_results_test()


# List of filenames and corresponding plot properties
files_and_properties = [
    {"filename": "Diamond_OPLS4.csv", "color": "red", "label": "Diamond_OPLS4", "marker": 'o', "markersize": 6},
    {"filename": "MCMM_OPLS4.csv", "color": "blue", "label": "MCMM_OPLS4", "marker": 's', "markersize": 8},
    {"filename": "MCLMCS_OPLS4.csv", "color": "yellow", "label": "MCLMCS_OPLS4", "marker": '^', "markersize": 8},
    {"filename": "SPMC_OPLS4.csv", "color": "purple", "label": "SPMC_OPLS4", "marker": 'D', "markersize": 8},
    {"filename": "LMCS_OPLS4.csv", "color": "green", "label": "LMCS_OPLS4", "marker": 'x', "markersize": 8}
]

# Initialize the merged dataframe with the first CSV file
initial_file = files_and_properties[0]["filename"]
merged_df = pd.read_csv(initial_file)
merged_df.dropna(subset=["Molecule", "Total_Running_Time"], inplace=True)
merged_df.rename(columns={"Total_Running_Time": f"Total_Running_Time_{files_and_properties[0]['label']}"}, inplace=True)

# Read the other CSV files and merge them into the merged dataframe
for file_prop in files_and_properties[1:]:
    filename = file_prop["filename"]
    df = pd.read_csv(filename)
    df.dropna(subset=["Molecule", "Total_Running_Time"], inplace=True)
    df.rename(columns={"Total_Running_Time": f"Total_Running_Time_{file_prop['label']}"}, inplace=True)
    merged_df = merged_df.merge(df, on="Molecule", how="inner")

# Plot the graph
plt.figure(figsize=(20, 8))

# Plot each line using the specified properties
for file_prop in files_and_properties:
    label = file_prop["label"]
    plt.plot(
        merged_df["Molecule"],
        merged_df[f"Total_Running_Time_{label}"],
        color=file_prop["color"],
        label=label,
        marker=file_prop["marker"],
        markersize=file_prop["markersize"],
        linestyle=':'
    )

plt.xlabel("Molecule")
plt.ylabel("Total Running Time(s) in logarithmic scale")
plt.title("Comparison of Algorithms")
plt.yscale("log")  # Set y-axis to logarithmic scale
plt.legend()
plt.grid(True)
plt.show()

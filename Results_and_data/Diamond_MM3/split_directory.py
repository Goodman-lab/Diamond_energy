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
import csv
import pandas as pd
from datetime import datetime

# The molecular name in acyclic alkane test set
# One can change the following list to target test set
name=["c4", "c5", "c5n", "c6", "c6n", "c7", "c7c", "c7n", "c8", "c8a",
      "c8s", "c8si", "c9", "c9n", "c9s", "c10n", "c12p", "c12p0", "c12p1", "c15",
      "c15helix", "c15a", "c15b", "c15c", "c15d", "c15e", "c15f", "c15g", "c16", "c17",
      "c17a", "c17l", "c18", "c18n", "c19", "c20n", "c21", "c26"]
total_number=len(name)

# The complete test molecule InChi information
#with open("Test_list.txt", "r") as f:
#    inchi_list=[]
#    for line in f.read().split():
#        print(line)
#        inchi_list.append(line)


#def count_files_in_directory(directory, extension):
#    return len([f for f in os.listdir(directory) if f.endswith(extension)])
# This function now returns a list of filenames with the given extension
def count_files_in_directory(directory, extension):
    return [f for f in os.listdir(directory) if f.endswith(extension)]


#dir_path='.'
#file_extension='.mae'
#file_count=count_files_in_directory(dir_path, file_extension)

# Generate 3D structure from InChi
# For molecules tested need to run maestro calculation,
# Here are the lines to turn them into .sdf format for MacroModel running input
def inchi_to_sdf():
    with open("Test_list.txt", "r") as f:
        #ßprint(f.read().split())
        inchi_list=[]
        for line in f.read().split():
            print(line)
            inchi_list.append(line)
        mol_totalnum=len(inchi_list)
        print("mollist_length", mol_totalnum)
        inchi_n=0
        for mol_inchi in inchi_list:
            mol=Chem.MolFromInchi(mol_inchi)
            rdkit.Chem.Draw.MolToFile(mol, name[inchi_n] + "_graph.png", size=(300, 300))
            mol_maestro=Chem.AddHs(mol, addCoords=True)
            AllChem.EmbedMolecule(mol_maestro)
            mol_maetsro_file=Chem.MolToMolFile(mol_maestro, name[inchi_n]+"_maestro.sdf")
            inchi_n+=1

# Generate input mol file for running MacroModel
# The function can use structconvert.py in Schrodinger/version/utilities/ directory
# To change .sdf format mol file to .mae MacroModel molecular structure format as input
# More instructions can be found in macromodel_reference_manual
def sdf2mae():
    i=0
    while i < total_number:
        # The below directory to call structconvert.py may vary in different computer
        command="/shared/shared/schrodinger/2022-1/utilities/structconvert "+name[i]+".sdf "+name[i]+"_maestro.mae"
        os.system(command)
        i+=1

# Generate command file for running MacroModel
# to run MacroModel from command line for batch running, instead of using its graphic interaction software Maestro by hand
# We need a molecular structure file in .mae format, which can be converted using structconvert.py in Schrodinger/version/utilities/ directory
# and a command file in .com format, which includes all the tasks (operation codes and instructions can be found in macromodel_reference_manual)
# function for generating MacroModel command .com file for corresponding molecular structure .mae file

# MCMM_OPLS4 template
def generate_comfile():
    num=0
    while num < total_number:
        with open(name[num]+"_MCMM_OPLS4.com", "w+") as f:
            f.write(name[num]+"_maestro.mae")
            f.write("\n")
            f.write(name[num]+"_MCMM_OPLS4-out.maegz")
            f.write("\n")
            f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
            f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
            f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
            f.write(" MCMM  100000      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" NANT       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" MCNV       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" MCSS       2      0      0      0    25.0000     0.0000     0.0000     0.0000\n")
            f.write(" MCOP       1      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" DEMX       0    166      0      0    25.0000    50.0000     0.0000     0.0000\n")
            f.write(" COMP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" MSYM       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" CHIG       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" AUOP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" CONV       2      0      0      0     0.0010     0.0000     0.0000     0.0000\n")
            f.write(" AUTO       0      3      1      1     1.0000     0.0000     0.0000     0.0000\n")
            f.write(" MINI       4      0   1000      0     0.0000     0.0000     0.0000     0.0000\n")
        f.close()
        num+=1
'''
# SPMC_MM3 template
def generate_comfile():
    num=0
    while num < total_number:
        with open(name[num]+"_SPMC_MM3.com", "w+") as f:
            f.write(name[num]+"_maestro.mae")
            f.write("\n")
            f.write(name[num]+"_SPMC_MM3-out.maegz")
            f.write("\n")
            f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" FFLD       2      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
            f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
            f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
            f.write(" SPMC  100000      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" NANT       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" MCNV       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" MCSS       2      0      0      0    25.0000     0.0000     0.0000     0.0000\n")
            f.write(" MCOP       1      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" DEMX       0    166      0      0    25.0000    50.0000     0.0000     0.0000\n")
            f.write(" COMP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" MSYM       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" CHIG       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" AUOP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" CONV       2      0      0      0     0.0010     0.0000     0.0000     0.0000\n")
            f.write(" AUTO       0      3      1      1     1.0000     0.0000     0.0000     0.0000\n")
            f.write(" MINI       4      0   1000      0     0.0000     0.0000     0.0000     0.0000\n")
        f.close()
        num+=1



# LMCS_OPLS4 template
def generate_comfile():
    num=0
    while num < total_number:
        with open(name[num]+"_LMCS_OPLS4.com", "w+") as f:
            f.write(name[num]+"_maestro.mae")
            f.write("\n")
            f.write(name[num]+"_LMCS_OPLS4-out.maegz")
            f.write("\n")
            f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
            f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
            f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
            f.write(" LMCS  100000      0      0      0     0.0000     0.0000     3.0000     6.0000\n")
            f.write(" NANT       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" MCSS       2      0      0      0    25.0000     0.0000     0.0000     0.0000\n")
            f.write(" MCOP       1      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" DEMX       0    166      0      0    25.0000    50.0000     0.0000     0.0000\n")
            f.write(" COMP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" MSYM       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" CHIG       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" AUOP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" CONV       2      0      0      0     0.0010     0.0000     0.0000     0.0000\n")
            f.write(" AUTO       0      3      1      1     1.0000     0.0000     0.0000     0.0000\n")
            f.write(" MINI       4      0   1000      0     0.0000     0.0000     0.0000     0.0000\n")
        f.close()
        num+=1
        
# MCLMCS_MM3 template
def generate_comfile():
    num=0
    while num < total_number:
        with open(name[num]+"_MCLMCS_MM3.com", "w+") as f:
            f.write(name[num]+"_maestro.mae")
            f.write("\n")
            f.write(name[num]+"_MCLMCS_MM3-out.maegz")
            f.write("\n")
            f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" FFLD       2      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
            f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
            f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
            f.write(" LMCS  100000      0      0      0     0.0000     0.0000     3.0000     6.0000\n")
            f.write(" NANT       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" MCNV       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" MCSS       2      0      0      0    25.0000     0.0000     0.0000     0.0000\n")
            f.write(" MCOP       1      0      0      0     0.5000     0.0000     0.0000     0.0000\n")
            f.write(" DEMX       0    166      0      0    25.0000    50.0000     0.0000     0.0000\n")
            f.write(" COMP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" MSYM       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" CHIG       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" AUOP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" CONV       2      0      0      0     0.0010     0.0000     0.0000     0.0000\n")
            f.write(" AUTO       0      3      1      1     1.0000     0.0000     0.0000     0.0000\n")
            f.write(" MINI       4      0   1000      0     0.0000     0.0000     0.0000     0.0000\n")
        f.close()
        num+=1
'''

# Submit batch MacroModel CS running
def run_CS_batch():
    a=0
    while a < total_number:
        command="bmin -WAIT "+name[a]+"_MCMM_OPLS4"
        os.system(command)
        a+=1

def run_diamond():
    # need to adjust for batch run later on
    command="python Diamond_energy.py 'InChI=1S/C18H38/c1-3-5-7-9-11-13-15-17-18-16-14-12-10-8-6-4-2/h3-18H2,1-2H3'"
    os.system(command)

def call_sdftomae():
    command="bash sdftomae.sh"
    os.system(command)

def generate_minifile4diamond():
    num=0
    while num < file_count:
        with open(str(num)+"_diamond_minimization.com", "w+") as f:
            f.write(str(num)+"_diamond.mae")
            f.write("\n")
            f.write(str(num)+"_diamond_minimization-out.maegz")
            f.write("\n")
            f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" DEBG      55      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
            f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
            f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
            f.write(" BGIN       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" CONV       2      0      0      0     0.0010     0.0000     0.0000     0.0000\n")
            f.write(" MINI       4      0   1000      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" END        0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
        f.close()
        num+=1

# Submit batch minimization running for Diamond_output
def run_mini_batch():
    b=0
    while b < total_number:
        command="bmin "+str(b)+"_diamond_minimization"
        os.system(command)
        b+=1

def split_directory():
    mol_order=30
    while mol_order < total_number:
        #command="mkdir "+name[mol_order]
        #os.system(command)
        #Best write it down and then can check the structure later on
        #with open("./"+name[mol_order]+"/"+name[mol_order]+".txt", "w+") as f:
        #    f.write(inchi_list[mol_order])
        #    f.write("\n")
        #f.close()
        command1="cp Diamond_energy.py ./"+name[mol_order]
        os.system(command1)
        #command2="cp generate_mol.py ./"+name[mol_order]
        #os.system(command2)
        #command3="cp sdftomae.sh ./"+name[mol_order]
        #os.system(command3)
        mol_order+=1


def run_diamondscript_parallel():
    mol_test=19
    while mol_test < total_number:
        # currently I have not tested long alkane chains, skip 17l, 18n ,20n
        if mol_test in [28, 29, 31, 32, 33, 34, 35, 36, 37]:
            mol_test+=1
            continue
        else:
            current_running_path="./"+name[mol_test]+"/"
            command_run="nohup python -u "+current_running_path+"generate_mol.py "+name[mol_test]+"> "+current_running_path+"running_record.log 2>&1 &"
            os.system(command_run)
        mol_test+=1


def read_results():
    for dir_name in os.listdir('.'):
        if os.path.isdir(dir_name):
            conf_num = None
            total_energy = None
            saddle_point = None
            multiple_minimum = None
            true_minimum = None
            lowest_energy = None
            saddle_point_count = 0
            multiple_minimum_count = 0
            true_minimum_count = 0

            directory = dir_name  # subdirectory to look
            extension = "diamondH_fretest.log"
            finalfretest_files = count_files_in_directory(directory, extension)

            csv_file_path = os.path.join(directory, f"{dir_name}.csv")
            with open(csv_file_path, 'w', newline='') as csvfile:
                csvwriter = csv.writer(csvfile)
                csvwriter.writerow(
                    ["Conf_num", "Total_Energy(kcal/mol)", "Saddle_point", "Multiple_minimum", "True_minimum",
                     "Lowest_energy", "Saddle_point_count", "Multiple_minimum_count", "True_minimum_count"])

                for file_name in finalfretest_files:
                    with open(os.path.join(directory, file_name), "r") as f:
                        energy_list = []
                        conf_num = file_name.split("_")[1]
                        for line in f:
                            if "Total Energy =" in line:
                                total_energy = float(line.split("=")[1].split("kcal")[0].strip())
                                energy_list.append(total_energy)
                            elif "True minimum" in line:
                                true_minimum = 1
                                true_minimum_count += 1
                            elif "Multiple minimum" in line:
                                multiple_minimum = 1
                                multiple_minimum_count += 1
                            elif "Saddle point" in line:
                                saddle_point = 1
                                saddle_point_count += 1

                        lowest_energy = min(energy_list, default=None)
                        csvwriter.writerow(
                            [conf_num, total_energy, saddle_point, multiple_minimum, true_minimum, lowest_energy,
                             saddle_point_count, multiple_minimum_count, true_minimum_count])


#if __name__ == "__main__":
 #   read_results()

def final_results():
    # Initialize an empty DataFrame to hold the final results
    final_df = pd.DataFrame(
        columns=["Subdirectory", "Lowest_Energy", "True_Min_Count", "Multiple_Min_Count", "Minimum_Found",
                 "Saddle_Point_Count"])

    # Loop through all subdirectories
    for subdir, _, _ in os.walk('.'):
        if subdir != '.':  # Skip the current directory
            csv_files = [f for f in os.listdir(subdir) if f.endswith('.csv')]
            for csv_file in csv_files:
                # Read each CSV file into a DataFrame
                df = pd.read_csv(os.path.join(subdir, csv_file))

                # Check if the columns exist in the DataFrame before accessing them
                lowest_energy = df['Lowest_energy'].min() if 'Lowest_energy' in df.columns else None
                true_min_count = df['True_minimum'].sum() if 'True_minimum' in df.columns else None
                multiple_min_count = df['Multiple_minimum'].sum() if 'Multiple_minimum' in df.columns else None
                saddle_point_count = df['Saddle_point'].sum() if 'Saddle_point' in df.columns else None

                minimum_found = true_min_count + multiple_min_count if true_min_count is not None and multiple_min_count is not None else None

                # Append to final DataFrame
                final_df = final_df.append({
                    "Subdirectory": subdir,
                    "Lowest_Energy": lowest_energy,
                    "True_Min_Count": true_min_count,
                    "Multiple_Min_Count": multiple_min_count,
                    "Minimum_Found": minimum_found,
                    "Saddle_Point_Count": saddle_point_count
                }, ignore_index=True)

    # Save the final DataFrame to a CSV file
    final_df.to_csv("final_results.csv", index=False)

    # Call the function

def adjust_resultsformat():
    # Load existing DataFrame from the CSV file
    df = pd.read_csv('final_results.csv')

    # Rename 'Subdirectory' column to 'Molecular Name'
    df.rename(columns={'Subdirectory': 'Molecular Name'}, inplace=True)

    # Remove './' from the 'Molecular Name' column
    df['Molecular Name'] = df['Molecular Name'].str.replace('./', '', regex=False)

    # Convert specified columns to integer type (except 'Lowest_Energy')
    columns_to_int = ['True_Min_Count', 'Multiple_Min_Count', 'Minimum_Found', 'Saddle_Point_Count']

    for col in columns_to_int:
        # Only convert if the column exists in the DataFrame
        if col in df.columns:
            df[col] = df[col].fillna(0).astype(int)  # Fill NaN values with 0 before converting to int

    # Save updated DataFrame back to CSV
    df.to_csv('final_results_updated.csv', index=False)


def read_minitime():
    final_df = pd.DataFrame(columns=['Subdirectory', 'Total_Running_Time'])

    # Function to calculate running time from log lines
    def calculate_time(start_line, end_line):
        start_time_str = start_line.split("Starting Time")[-1].strip()
        end_time_str = end_line.split("termination")[-1].strip()

        start_time = datetime.strptime(start_time_str, '%d-%b-%Y %H:%M:%S')
        end_time = datetime.strptime(end_time_str, '%d-%b-%Y %H:%M:%S')

        return (end_time - start_time).seconds

    # Navigate through each subdirectory
    for dirpath, dirnames, filenames in os.walk('.'):
        total_time = 0

        # Check each file in the current subdirectory
        for filename in filenames:
            if filename.endswith("minimization.log"):
                with open(os.path.join(dirpath, filename), "r") as f:  # Fixed the undefined variables here
                    lines = f.readlines()
                    start_line_candidates = [line for line in lines if "Starting Time" in line]
                    end_line_candidates = [line for line in lines if "normal termination" in line]

                    if len(start_line_candidates) == 0 or len(end_line_candidates) == 0:
                        print(f"Skipping {filename} as it doesn't contain either Starting Time or Ending Time.")
                        continue

                    start_line = start_line_candidates[0]
                    end_line = end_line_candidates[0]
                    # Calculate time difference goes here
                    total_time += calculate_time(start_line, end_line)

        # Add data to the final DataFrame
        final_df = final_df.append({'Subdirectory': dirpath, 'Total_Running_Time': total_time}, ignore_index=True)

    # Save the DataFrame as a CSV file
    final_df.to_csv('total_running_times.csv', index=False)

#read_minitime()
adjust_resultsformat()
#final_results()

#inchi_to_sdf()
#generate_comfile()
#run_CS_batch()
#run_diamond()
#call_sdftomae()
#generate_minifile4diamond()
#run_mini_batch()
#split_directory()
#run_diamondscript_parallel()

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
import sys
import glob

def count_files_in_directory(directory, pattern):
    return len(glob.glob(os.path.join(directory, pattern)))

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


#def count_files_in_directory(directory, extension):
#    return len([f for f in os.listdir(directory) if f.endswith(extension)])
#dir_path=os.getcwd()
#file_extension='.sdf'
#file_count=count_files_in_directory(dir_path, file_extension)
#print(file_count)
def count_files_in_directory(directory, pattern):
    return len(glob.glob(os.path.join(directory, pattern)))

def files_in_directory(directory, pattern):
    return [os.path.basename(f) for f in glob.glob(os.path.join(directory, pattern))]




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
            #AllChem.EmbedMolecule(mol_maestro)
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

def generate_comfile(method_name):
    if method_name=="MCMM_OPLS4":
        # MCMM_OPLS4 template
        num=0
        while num < total_number:
            with open(script_directory+name[num]+"_MCMM_OPLS4.com", "w+") as f:
                f.write(script_directory+name[num]+"_maestro.mae")
                f.write("\n")
                f.write(script_directory+name[num]+"_MCMM_OPLS4-out.mae")
                f.write("\n")
                f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
                f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
                f.write(" MCMM   10000      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
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
                f.write(" MINI       1      0   5000      0     0.0000     0.0000     0.0000     0.0000\n")
            f.close()
            command="bmin -WAIT {0}{1}_MCMM_OPLS4".format(script_directory, name[num])
            os.system(command)
            num+=1
    elif method_name=="MCMM_MM3":
        # MCMM_MM3 template
        num=0
        while num < total_number:
            with open(script_directory+name[num]+"_MCMM_MM3.com", "w+") as f:
                f.write(script_directory+name[num]+"_maestro.mae")
                f.write("\n")
                f.write(script_directory+name[num]+"_MCMM_MM3-out.mae")
                f.write("\n")
                f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" FFLD       2      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
                f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
                f.write(" MCMM   10000      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
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
                f.write(" MINI       1      0   5000      0     0.0000     0.0000     0.0000     0.0000\n")
            f.close()
            command="bmin -WAIT {0}{1}_MCMM_MM3".format(script_directory, name[num])
            os.system(command)
            num+=1
    elif method_name=="SPMC_MM3":
        # SPMC_MM3 template
        num=0
        while num < total_number:
            with open(script_directory+name[num]+"_SPMC_MM3.com", "w+") as f:
                f.write(script_directory+name[num]+"_maestro.mae")
                f.write("\n")
                f.write(script_directory+name[num]+"_SPMC_MM3-out.mae")
                f.write("\n")
                f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" FFLD       2      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
                f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
                f.write(" SPMC   10000      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
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
                f.write(" MINI       1      0   5000      0     0.0000     0.0000     0.0000     0.0000\n")
            f.close()
            command="bmin -WAIT {0}{1}_SPMC_MM3".format(script_directory, name[num])
            os.system(command)
            num+=1
    elif method_name=="SPMC_OPLS4":
        # SPMC_OPLS4 template
        num=0
        while num < total_number:
            with open(script_directory+name[num]+"_SPMC_OPLS4.com", "w+") as f:
                f.write(script_directory+name[num]+"_maestro.mae")
                f.write("\n")
                f.write(script_directory+name[num]+"_SPMC_OPLS4-out.mae")
                f.write("\n")
                f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
                f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
                f.write(" SPMC   10000      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
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
                f.write(" MINI       1      0   5000      0     0.0000     0.0000     0.0000     0.0000\n")
            f.close()
            command="bmin -WAIT {0}{1}_SPMC_OPLS4".format(script_directory, name[num])
            os.system(command)
            num+=1
    elif method_name=="LMCS_OPLS4":
        # LMCS_OPLS4 template
        num=0
        while num < total_number:
            with open(script_directory+name[num]+"_LMCS_OPLS4.com", "w+") as f:
                f.write(script_directory+name[num]+"_maestro.mae")
                f.write("\n")
                f.write(script_directory+name[num]+"_LMCS_OPLS4-out.mae")
                f.write("\n")
                f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
                f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
                f.write(" LMCS   10000      0      0      0     0.0000     0.0000     3.0000     6.0000\n")
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
                f.write(" MINI       1      0   5000      0     0.0000     0.0000     0.0000     0.0000\n")
            f.close()
            command="bmin -WAIT {0}{1}_LMCS_OPLS4".format(script_directory, name[num])
            os.system(command)
            num+=1
    elif method_name=="LMCS_MM3":
        # LMCS_MM3 template
        num=0
        while num < total_number:
            with open(script_directory+name[num]+"_LMCS_MM3.com", "w+") as f:
                f.write(script_directory+name[num]+"_maestro.mae")
                f.write("\n")
                f.write(script_directory+name[num]+"_LMCS_MM3-out.mae")
                f.write("\n")
                f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" FFLD       2      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
                f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
                f.write(" LMCS   10000      0      0      0     0.0000     0.0000     3.0000     6.0000\n")
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
                f.write(" MINI       1      0   5000      0     0.0000     0.0000     0.0000     0.0000\n")
            f.close()
            command="bmin -WAIT {0}{1}_LMCS_MM3".format(script_directory, name[num])
            os.system(command)
            num+=1
    elif method_name=="MCLMCS_MM3":
        # MCLMCS_MM3 template
        num=0
        while num < total_number:
            with open(script_directory+name[num]+"_MCLMCS_MM3.com", "w+") as f:
                f.write(script_directory+name[num]+"_maestro.mae")
                f.write("\n")
                f.write(script_directory+name[num]+"_MCLMCS_MM3-out.mae")
                f.write("\n")
                f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" FFLD       2      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
                f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
                f.write(" LMCS   10000      0      0      0     0.0000     0.0000     3.0000     6.0000\n")
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
                f.write(" MINI       1      0   5000      0     0.0000     0.0000     0.0000     0.0000\n")
            f.close()
            command="bmin -WAIT {0}{1}_MCLMCS_MM3".format(script_directory, name[num])
            os.system(command)
            num+=1
    elif method_name=="MCLMCS_OPLS4":
        # MCLMCS_OPLS4 template
        num=0
        while num < total_number:
            with open(script_directory+name[num]+"_MCLMCS_OPLS4.com", "w+") as f:
                f.write(script_directory+name[num]+"_maestro.mae")
                f.write("\n")
                f.write(script_directory+name[num]+"_MCLMCS_OPLS4-out.mae")
                f.write("\n")
                f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
                f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
                f.write(" LMCS   10000      0      0      0     0.0000     0.0000     3.0000     6.0000\n")
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
                f.write(" MINI       1      0   5000      0     0.0000     0.0000     0.0000     0.0000\n")
            f.close()
            command="bmin -WAIT {0}{1}_MCLMCS_OPLS4".format(script_directory, name[num])
            os.system(command)
            num+=1


def CSoutput_split2singlefile():
    command="bash {}results_split.sh".format(script_directory)
    os.system(command)

def frequency_test():
    num=0
    while num < file_count:
        with open(script_directory+minifiles_name[num]+"_fretest.com", "w+") as f:
            f.write(script_directory+minifiles_name[num]+".mae")
            f.write("\n")
            f.write(script_directory+minifiles_name[num]+"_fretest-out.mae")
            f.write("\n")
            f.write(" DEBG     211      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" FFLD       2      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
            f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
            f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" MINI       4      0   5000      0     0.0000     0.0000     0.1000     0.0000\n")
            f.write(" ELST       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" MTST       0      0      0      0     2.0000     0.0000     0.0000     0.0000\n")
            f.write(" RRHO       0      0      1      0   300.0000     1.0000     2.0000     0.0000\n")
        f.close()
        command="bmin -WAIT {0}{1}_fretest".format(script_directory, minifiles_name[num])
        os.system(command)
        num+=1



script_path=os.path.join(os.getcwd(), sys.argv[1])
script_directory=script_path+"/"
dir_path=script_directory
method=sys.argv[1]
generate_comfile(method)
CSoutput_split2singlefile()
minimolfile_pattern='*-mini-*.mae'# Your file pattern
file_count=count_files_in_directory(dir_path, minimolfile_pattern)
print("file_count:", file_count)
file_outputmae=files_in_directory(dir_path, minimolfile_pattern)
minifiles_name=[]
for file in file_outputmae:
    minifiles_name.append(file.split(".")[0])
minifiles_name.sort()
print(minifiles_name)
frequency_test()


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

# The molecular name in extension work test set
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


#dir_path='.'
#file_extension='.mae'
#file_count=count_files_in_directory(dir_path, file_extension)


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


def split_directory():
    test_order=0
    while test_order < total_method:
        command="mkdir "+Maestro_test[test_order]
        os.system(command)
        command1="cp maestro_mol.py ./"+Maestro_test[test_order]
        os.system(command1)
        command2="cp results_split.sh  ./"+Maestro_test[test_order]
        os.system(command2)
        command3="cp *.mae  ./" + Maestro_test[test_order]
        os.system(command3)
        test_order+=1

def run_maestroscript_parallel():
    test_order=0
    while test_order < total_method:
        current_running_path="./"+Maestro_test[test_order]+"/"
        command_run="nohup python -u "+current_running_path+"maestro_mol.py "+Maestro_test[test_order]+"> "+current_running_path+"running_record_2.log 2>&1 &"
        os.system(command_run)
        test_order+=1





#split_directory()
run_maestroscript_parallel()

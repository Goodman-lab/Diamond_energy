# Diamond Energy Algorithm Experiments

This repository contains scripts and files for running experiments using the Diamond Energy algorithm in batch mode.

## Contents

1. **Test_list.txt**
- Contains the InChI representations of test molecules used in these experiments.

2. **sdftomae.sh**
- Converts file formats during the process.
- Generates conformer output files (.sdf) from the Diamond Energy algorithm.
- These conformer files can then be used in the MacroModel suite for minimization (requires .mae format).

3. **split_directory.py**
- The initial script to run.
- Generates subdirectories for each test molecule.
- Each subdirectory is named after the corresponding test molecule.
- Automatically copies the necessary files into each subdirectory.
- Initiates the conformer search process for each test molecule within its designated subdirectory.

4. **generate_mol.py**
- Calls the diamondenergy.py script for each test molecule in its respective subdirectory.
- Controls the conformer search process within that subdirectory.


5. **diamondenergy.py**
- The core script implementing the Diamond Energy algorithm.

6. **read_energyvalue.py**
- After completing all processes, this script reads energy value results from all subdirectories.
- Stores the results in a single .csv file in the main directory.

7. **read_runningtime.py**
- Similar to the previous script, this one reads running time results from all subdirectories.
- Also stores the results in a single .csv file in the main directory.

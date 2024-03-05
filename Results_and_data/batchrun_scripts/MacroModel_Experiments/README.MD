# Conformer Search Algorithms in MacroModel Experiments

This repository contains scripts and files for running four leading conformer search algorithms (**SPMC**, **LMCS**, **MCLMCS**, and **MCMM**) in batch mode using the MacroModel suite. Please note that these scripts require a valid Schrödinger license to execute.

## Contents

1. **Test_list.txt**
    - Contains the InChI representations of test molecules used in these experiments.

2. **Results_split.sh**
    - Splits the compact minimization results for a specific test molecule (ending in `.maegz`) into individual files (ending in `.mae`).
    - These split files can then be used for frequency tests conducted in MacroModel.

3. **Split_directory_mae.py**
    - The initial script to run.
    - Generates subdirectories for each test algorithm.
    - Each subdirectory is named after the corresponding test algorithm and the force field used.
    - Automatically copies the required files into each subdirectory.
    - Initiates the conformer search process for all test molecules within their designated subdirectories, using the corresponding conformer search method and force field.

4. **maestro_mol.py**
    - Executes each conformer search method and force field for all test molecules in their respective subdirectories.
    - Controls the conformer search process within each subdirectory.

5. **Read_energyvalue.py**
    - After completing all processes, this script reads energy value results from all subdirectories.
    - Stores the results in a single `.csv` file in the main directory.

6. **Read_runningtime.py**
    - Similar to the previous script, this one reads running time results from all subdirectories.
    - Also stores the results in a single `.csv` file in the main directory.


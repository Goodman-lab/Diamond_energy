# *Diamond Energy*

Jonathan M. Goodman, Mengman Wei

A systematic conformational searching algorithm on a diamond lattice framework

#################################################################

# CONTENTS

1. Release Notes

2. Requirements and Setup

3. Usage


#################################################################

# Release Notes

*Diamond Energy* is an algorithm developed in Python 3.6+ and RDKit environment for conducting systematic conformational searches of molecules. 

The algorithm operates in a series of steps:

Input Parsing: Initially, the algorithm takes the InChI of the target molecule as input to extract its molecular information.

Framework Assumption: The algorithm uses the premise that low-energy conformations of a molecule are likely to closely resemble diamond lattice structures. This assumption enables a more efficient search.

Conformer Generation: Utilizing tree-search algorithm, *Diamond Energy* generates potential low-energy conformers of the molecule. These conformers are designed to fit within the constraints of a diamond lattice framework.

Refinement: Finally, the output generated conformers could be further minimized using Quantum Mechanics/Molecular Mechanics (QM/MM) methods.

This methodical approach allows *Diamond Energy* to rapidly identify all potential low-energy conformers of a given molecule, aiding in subsequent computational or experimental analysis.

#################################################################

# Requirements and Setup

The script is compatible with Python 3.6 or higher and the RDKit 2020.09.1

Development and testing were conducted on Linux and MacOS operating systems.

#################################################################

# Getting Started

Make sure Python 3.6+ is installed on your machine.
RDKit should be installed and properly configured. If not, installation instructions can be found on the RDKit website(https://www.rdkit.org/docs/Install.html).

1. Download the Diamond_energy.py script.
2. Navigate to the directory where the script is located using the terminal. 
3. Activate RDKit enviroment and run the program from the terminal.

# Correct Usage Syntax

__Global Search Mode__

To perform a comprehensive conformational search:

python Diamond_energy.py '<Molecular_InChI>' 

__Specifying the Number of Conformers__

To generate and display a specific number of conformers:

python Diamond_energy.py '<Molecular_InChI>'  '<Number_of_Conformers>'

# Example: Heptane

To perform a conformational search for heptane, use the following command:

python Diamond_energy.py "InChI=1S/C7H16/c1-3-5-7-6-4-2/h3-7H2,1-2H3"

To obtain the first conformer of the heptane, use the following command:

python Diamond_energy.py "InChI=1S/C7H16/c1-3-5-7-6-4-2/h3-7H2,1-2H3" 0

# Output

The program automatically conducts a conformational search for the specified molecule, such as heptane in the provided example. The terminal will present various details about the search process and its results. The output may include, but is not limited to, the following:

Lowest Energy Conformation: The conformation with the minimum energy.

Energies of Each Conformation: Energy values associated with each identified conformation.

Total Number of Accessible Conformations: Count of conformations that have accessible energy levels.

...and more.

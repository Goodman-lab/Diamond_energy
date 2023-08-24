# Diamond_energy

Jonathan M. Goodman, Mengman Wei

A systematic conformational searching algorithm on a diamond lattice framework

#################################################################

# CONTENTS

1. Release Notes

2. Requirements and Setup

3. Usage


#################################################################

# Release Notes

Diamond_energy is a systematic conformational searching algorithm written in Python 3.6+. It operates on the premise that local minima in the energy landscape closely resemble diamond lattice structures. By leveraging this assumption, the program rapidly identifies all low-energy conformations through adjustments to the molecule's internal coordinates within the constraints of a diamond lattice framework.

The program initiates by extracting molecular information from the provided InChI input. It then generates conformers of the molecule within a diamond lattice framework, employing internal coordinate tree-search algorithms. These low-energy conformers output, optimized to fit the diamond lattice structure, serve as starting points for subsequent minimization via Quantum Mechanics/Molecular Mechanics (QM/MM) methods.

#################################################################

# Requirements and Setup

The script is compatible with Python 3.6 or higher and requires the RDKit library.

Development and testing were conducted on Linux and MacOS operating systems.

#################################################################

# Getting Started

1. Download the Diamond_energy.py script.
2. Navigate to the directory where the script is located using the terminal.
3. Run the program from the terminal.

# Correct Usage Syntax

__Global Search Mode__

To perform a comprehensive conformational search:

python Diamond_energy.py <MolecularInChI>

__Specifying the Number of Conformers__

To generate and display a specific number of conformers:

python Diamond_energy.py <MolecularInChI> <Number_of_Conformers>

# Example: Heptane

To perform a conformational search for heptane, use the following command:

python Diamond_energy.py "InChI=1S/C7H16/c1-3-5-7-6-4-2/h3-7H2,1-2H3"

# Output

The program automatically conducts a conformational search for the specified molecule, such as heptane in the provided example. The terminal will present various details about the search process and its results. The output may include, but is not limited to, the following:

Lowest Energy Conformation: The conformation with the minimum energy.

Energies of Each Conformation: Energy values associated with each identified conformation.

Total Number of Accessible Conformations: Count of conformations that have accessible energy levels.

Low-Energy Conformers 3D File: Conformers with low energies are saved in a .sdf file.

...and more.

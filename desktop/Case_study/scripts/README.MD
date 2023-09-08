# Scripts Introduction

This repository contains a set of Python and shell scripts used for conformer generation and analysis in case study. Below is a breakdown of each script's functionality.

---

## Table of Contents

1. [Head_control.py](#head_controlpy)
2. [Diamond_energy_main_algorithm.py](#diamond_energy_main_algorithmpy)
3. [generate_mol.py](#generate_molpy)
4. [Shell Files (.sh)](#shell-files-sh)
5. [Count_diamondoutput_helix.py](#count_diamondoutput_helixpy)
6. [Count_diamondmini_helix.py](#count_diamondmini_helixpy)
7. [Data Analysis Scripts](#data-analysis-scripts)

---

## Head_control.py

### Purpose

This is the main driver script that automates the conformer search calculations for all molecules in the `Test_list.txt` file.

#### Functions

- `split_directory`: Creates subdirectories for each molecule listed in `Test_list.txt` and copies essential scripts into these directories.

##### Files Copied

- `Diamond_energy_main_algorithm.py`: Conformer generation
- `generate_mol.py`: Prepares output for MacroModel minimization
- Shell files for format conversion (`maetosdf.sh`, `sdftomae.sh`)

---

## Diamond_energy_main_algorithm.py

### Purpose

This is the core script containing the algorithm for conformer generation based on the diamond lattice framework.

---

## generate_mol.py

### Purpose

Processes the output from `Diamond_energy_main_algorithm.py`.

#### Functions

- Adds hydrogens to the structure
- Converts `.sdf` format to a MacroModel-readable format
- Submits minimization tasks to MacroModel

> Note: Users can change the force field, minimization algorithm, and other parameters. (Refer to the MacroModel manual for details.)

---

## Shell Files (.sh)

### Purpose

These shell scripts are used for converting molecular structure formats.

---

## Count_diamondoutput_helix.py

### Purpose

Reads `runrecord.log` from `Diamond_energy_main_algorithm.py` and calculates the helix ratio based on its output.

---

## Count_diamondmini_helix.py

### Purpose

Processes minimization output files (`.mae`) from MacroModel.

#### Functions

- Converts `.mae` to `.sdf` format
- Utilizes RDKit to calculate dihedral angles
- Computes the helix ratio

---

## Data Analysis Scripts

### Files

- `Data_analysis.py`
- `Data_analysis(R/L).py`

### Purpose

Generate dot plots for visualizing output data.


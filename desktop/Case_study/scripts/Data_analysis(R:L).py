import pandas as pd
import matplotlib.pyplot as plt

# Function to convert string ratio to float
def ratio_to_float(ratio_str):
    x, y = map(float, ratio_str.split(":"))
    return x / y

# Read data from the three files into dataframes
molecule_df = pd.read_csv("molecule_helix_ratios_nosolvent.csv")
helix_df = pd.read_csv("helix_ratios.csv")
solvent_df = pd.read_csv("molecule_helix_ratios_solvent.csv")

# Merge the three dataframes based on "Molecular Name" column
merged_df = pd.merge(molecule_df, helix_df, on="Molecular Name", suffixes=("_molecule", "_helix"))
merged_df = pd.merge(merged_df, solvent_df, on="Molecular Name", suffixes=("", "_solvent"))

# Convert the "Formatted Helix Ratio" columns to a float
merged_df["Float Ratio_molecule"] = merged_df["Formatted Helix Ratio_molecule"].apply(ratio_to_float)
merged_df["Float Ratio_helix"] = merged_df["Formatted Helix Ratio_helix"].apply(ratio_to_float)
merged_df["Float Ratio_solvent"] = merged_df["Formatted Helix Ratio"].apply(ratio_to_float)

# Sort the merged dataframe
merged_df = merged_df.sort_values(by="Float Ratio_molecule")

# Extract relevant columns for plotting
molecule_names = merged_df["Molecular Name"]
molecule_float_ratios = merged_df["Float Ratio_molecule"]
helix_float_ratios = merged_df["Float Ratio_helix"]
solvent_float_ratios = merged_df["Float Ratio_solvent"]
molecule_ratios = merged_df["Formatted Helix Ratio_molecule"]
helix_ratios = merged_df["Formatted Helix Ratio_helix"]
solvent_ratios = merged_df["Formatted Helix Ratio"]

# Create figure and axis
fig, ax = plt.subplots()

# Plot data
ax.plot(molecule_names, molecule_float_ratios, marker='o', linestyle='-', label="DiamondEnergy_MM3(gas)")
ax.plot(molecule_names, helix_float_ratios, marker='s', linestyle='--', label="DiamondEnergy_output")
ax.plot(molecule_names, solvent_float_ratios, marker='^', linestyle='-.', label="DiamondEnergy_MMFFs(CHCl3)")

# Determine y-axis range
min_y = min(min(molecule_float_ratios), min(helix_float_ratios), min(solvent_float_ratios))
max_y = max(max(molecule_float_ratios), max(helix_float_ratios), max(solvent_float_ratios))
margin = (max_y - min_y) * 0.1

# Set y-axis limits
ax.set_ylim(min_y - margin, max_y + margin)

# Annotate each point
for i, (name, ratio) in enumerate(zip(molecule_names, molecule_ratios)):
    ax.annotate(ratio, (name, molecule_float_ratios.iloc[i]), textcoords="offset points", xytext=(0,5), ha='center')

for i, (name, ratio) in enumerate(zip(molecule_names, helix_ratios)):
    ax.annotate(ratio, (name, helix_float_ratios.iloc[i]), textcoords="offset points", xytext=(0,-25), ha='center')

for i, (name, ratio) in enumerate(zip(molecule_names, solvent_ratios)):
    ax.annotate(ratio, (name, solvent_float_ratios.iloc[i]), textcoords="offset points", xytext=(0,-10), ha='center')

# Rotate x-axis labels
plt.xticks(rotation=90)

# Add labels and legend
ax.set_xlabel("Molecular Name")
ax.set_ylabel("Float Helix Ratio")
ax.set_title("Right helix:Left helix results")
ax.legend()

# Show plot
plt.tight_layout()
plt.show()



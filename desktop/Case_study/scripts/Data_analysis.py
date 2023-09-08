import pandas as pd
import matplotlib.pyplot as plt

# Read data from "molecule_helix_ratios_nosolvent.csv", "helix_ratios.csv", and "molecule_helix_ratios_solvent.csv" into dataframes
molecule_df = pd.read_csv("molecule_helix_ratios_nosolvent.csv")
helix_df = pd.read_csv("helix_ratios.csv")
solvent_df = pd.read_csv("molecule_helix_ratios_solvent.csv")

# Merge the dataframes based on "Molecular Name" column
merged_df = pd.merge(molecule_df, helix_df, on="Molecular Name", suffixes=("_molecule", "_helix"))
merged_df = pd.merge(merged_df, solvent_df, on="Molecular Name", suffixes=("", "_solvent"))

# Print the columns in merged_df to debug
print("Columns in merged_df: ", merged_df.columns)

# Extract "Total Helix" columns
molecule_total_helix = merged_df["Total Helix Ratio_molecule"]
helix_total_helix = merged_df["Total Helix Ratio_helix"]
solvent_total_helix = merged_df["Total Helix Ratio"]

# Sort the dataframe by "Total Helix Ratio_molecule" in ascending order
merged_df = merged_df.sort_values(by="Total Helix Ratio_molecule")

# Extract the relevant columns for plotting
molecule_names = merged_df["Molecular Name"]
molecule_ratios = merged_df["Formatted Helix Ratio_molecule"]
helix_ratios = merged_df["Formatted Helix Ratio_helix"]
solvent_ratios = merged_df["Formatted Helix Ratio"]

# Create a figure and axis
fig, ax = plt.subplots()

# Plot data without shifts (as you wanted the original data represented)
ax.plot(molecule_names, molecule_total_helix, marker='o', linestyle='-', linewidth=2, label="DiamondEnergy_MM3(gas)", zorder=3)
ax.plot(molecule_names, helix_total_helix, marker='s', linestyle='--', linewidth=2, label="DiamondEnergy_output", zorder=2)
ax.plot(molecule_names, solvent_total_helix, marker='^', linestyle='-.', linewidth=2, label="DiamondEnergy_MMFFs(CHCl3)", zorder=1)

# Annotate each point with the actual ratio in string format
for i, (name, ratio) in enumerate(zip(molecule_names, molecule_ratios)):
    ax.annotate(ratio, (name, molecule_total_helix.iloc[i]), textcoords="offset points", xytext=(0,10), ha='center', zorder=4)

for i, (name, ratio) in enumerate(zip(molecule_names, helix_ratios)):
    ax.annotate(ratio, (name, helix_total_helix.iloc[i]), textcoords="offset points", xytext=(0,-10), ha='center', zorder=4)

for i, (name, ratio) in enumerate(zip(molecule_names, solvent_ratios)):
    ax.annotate(ratio, (name, solvent_total_helix.iloc[i]), textcoords="offset points", xytext=(0, -10), ha='center', zorder=4)

# Rotate x-axis labels for better visibility
plt.xticks(rotation=90)

# Add labels and legend
ax.set_xlabel("Molecular Name")
ax.set_ylabel("Total Helix Ratio")
ax.set_title("Total Helix Ratio in Conformer Search Results")
ax.legend()

# Show the plot
plt.tight_layout()
plt.show()



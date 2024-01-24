import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Load the merged data
merged_df = pd.read_csv('merged_data.csv')

# Plotting with logarithmic scale and larger text
plt.figure(figsize=(14, 7))

# Increase font size globally for the plot
plt.rcParams.update({'font.size': 12})

# Diamond data in red
plt.plot(merged_df['Molecule'], merged_df['Diamond'], 'r--o', label='Diamond')

# SPMC data
plt.plot(merged_df['Molecule'], merged_df['SPMC'], 'b--o', label='SPMC')

# Other columns
for column in merged_df.columns[2:]:
    if column not in ['Diamond', 'SPMC']:
        plt.plot(merged_df['Molecule'], merged_df[column], '--o', label=column)

# Set the y-axis to a logarithmic scale to spread out the data points
plt.yscale('log')

# Customizing the plot with larger text
plt.xlabel('Molecule', fontsize=14)
plt.ylabel('Average Conformer Generate Time (log scale)', fontsize=14)
plt.title('Average Conformer Generate Time per Molecule (Log Scale)', fontsize=16)
plt.legend(fontsize=12)
plt.xticks(rotation=90, fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()

# Show the plot
plt.show()
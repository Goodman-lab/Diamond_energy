import re
import os
import pandas as pd
import glob

name=["c4", "c5", "c5n", "c6", "c6n", "c7", "c7c", "c7n", "c8", "c8a",
      "c8s", "c8si", "c9", "c9n", "c9s", "c10n", "c12p", "c12p0", "c12p1", "c15",
      "c15helix", "c15a", "c15b", "c15c", "c15d", "c15e", "c15f", "c15g", "c16", "c17",
      "c17a", "c17l", "c18", "c18n", "c19", "c20n", "c21", "c26"]

def read_eachMaestrotest():
    for test_method in Maestro_test:
        summary_df = pd.DataFrame(
            columns=['Molecule', 'Lowest_Total_Energy', 'True_minimum_count', 'Multiple_minimum_count',
                     'Saddle_point_count'])

        for molecule in name:
            energy_df = pd.DataFrame(
                columns=['Random_Order', 'Total_Energy', 'True_minimum', 'Multiple_minimum', 'Saddle_point'])

            log_files = glob.glob(f"./{test_method}/{molecule}_{test_method}-out-mini-*_fretest.log")
            valid_filename_pattern = re.compile(rf'^{molecule}_{test_method}-out-mini(-\d+)?_fretest\.log$')

            for log_file in log_files:
                filename = os.path.basename(log_file)

                if valid_filename_pattern.match(filename):
                    try:
                        if '-out-mini-' in filename:
                            random_order = int(filename.split('-out-mini-')[1].split('_fretest.log')[0])
                        else:
                            random_order = 1


                        with open(log_file, 'r', encoding='utf-8') as f:
                            lines = f.readlines()
                            total_energy = None
                            true_minimum = 0
                            multiple_minimum = 0
                            saddle_point = 0

                            for line in lines:
                                if "Total Energy" in line:
                                    total_energy = float(line.split('=')[-1].strip().split()[0])
                                elif "True minimum" in line:
                                    true_minimum = 1
                                elif "Multiple minimum" in line:
                                    multiple_minimum = 1
                                elif "Saddle point" in line:
                                    saddle_point = 1

                            new_row = pd.DataFrame({'Random_Order': [random_order], 'Total_Energy': [total_energy],
                                                    'True_minimum': [true_minimum],
                                                    'Multiple_minimum': [multiple_minimum],
                                                    'Saddle_point': [saddle_point]})
                            energy_df = pd.concat([energy_df, new_row], ignore_index=True)

                    except FileNotFoundError:
                        print(f"File {log_file} not found.")
                        continue
                else:
                    print(f"Skipping invalid filename: {filename}")

            energy_df.to_csv(os.path.join(".", test_method, f'{molecule}_fre.csv'), index=False)


            lowest_energy = energy_df['Total_Energy'].min()

            # Count conformers within specified energy ranges
            energy_within_1kcal = energy_df[energy_df['Total_Energy'] <= lowest_energy + 1].shape[0]
            energy_within_3kcal = energy_df[energy_df['Total_Energy'] <= lowest_energy + 3].shape[0]
            energy_within_5kcal = energy_df[energy_df['Total_Energy'] <= lowest_energy + 5].shape[0]

            true_minimum_count = energy_df['True_minimum'].sum()
            multiple_minimum_count = energy_df['Multiple_minimum'].sum()
            saddle_point_count = energy_df['Saddle_point'].sum()

            new_row = pd.DataFrame({'Molecule': [molecule], 'Lowest_Total_Energy': [lowest_energy],
                                    'Conformers_within_1kcal': [energy_within_1kcal],
                                    'Conformers_within_3kcal': [energy_within_3kcal],
                                    'Conformers_within_5kcal': [energy_within_5kcal],
                                    'True_minimum_count': [true_minimum_count],
                                    'Multiple_minimum_count': [multiple_minimum_count],
                                    'Saddle_point_count': [saddle_point_count]})
            summary_df = pd.concat([summary_df, new_row], ignore_index=True)
        summary_df.to_csv(os.path.join(".", test_method, 'totalfre.csv'), index=False)


# The comparing experiment names of Macromodel with algorithm and force field
Maestro_test = ["MCMM_OPLS4", "MCMM_MM3", "SPMC_OPLS4", "SPMC_MM3",
                "LMCS_OPLS4", "LMCS_MM3", "MCLMCS_OPLS4", "MCLMCS_MM3"]

# Execute the function
read_eachMaestrotest()

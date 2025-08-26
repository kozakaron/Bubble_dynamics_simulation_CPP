"""
This script is used to test the stability of the solver with different settings. 
It runs a large number of simulations with various mechanisms and parameters designed to give a hard time to the solver. 
There are a total of 5 bruteforce parameter studies, statistics are printted on the console and saved into a txt file.
Main indicator: rate of succesfull simulations [%] combined.

Usage: python3 ./dev/stability_test.py <txtname> <detailed_description>
       <txtname> is the name of the txt file where the results will be saved.
       <detailed_description> is a detailed description of the test.

E.g.: python3 ./dev/stability_test.py test_name "Here you can detail the settings of the solver."
      The results are saved in ./_parameter_studies/stability_test/test_name.txt
"""


import sys
import os
import shutil
import time
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from interface import python_interface as api

base_dir = './_parameter_studies/stability_test/'
semi_success_time = 1000e-6
timeout = 30.0

parameter_studies = dict(
    chemkin_ar_he__high_p_A = {
        'mechanism': 'chemkin_ar_he',
        'R_E': {'type': 'LinearRange', 'start': 0.5e-06, 'end': 20e-6, 'num_steps': 20},
        'ratio': {'type': 'Const', 'value': 1.0},
        'species': ['O2'],
        'fractions': [1.0],
        'P_amb': {'type': 'Const', 'value': 0.5e5},
        'T_inf': {'type': 'Const', 'value': 293.15},
        'alfa_M': {'type': 'Const', 'value': 0.35},
        'P_v': {'type': 'Const', 'value': 2338.1},
        'mu_L': {'type': 'Const', 'value': 0.001},
        'rho_L': {'type': 'Const', 'value': 998.2},
        'c_L': {'type': 'Const', 'value': 483.0},
        'surfactant': {'type': 'Const', 'value': 1.0},
        'enable_heat_transfer': True,
        'enable_evaporation': True,
        'enable_reactions': True,
        'enable_dissipated_energy': True,
        'target_specie': 'H2',
        'excitation_params': [
            {'type': 'LinearRange', 'start': -12e5, 'end': -22e5, 'num_steps': 5},
            {'type': 'LinearRange', 'start': 2500.0, 'end': 12500.0, 'num_steps': 5},
            {'type': 'Const', 'value': 1.0}
        ],
        'excitation_type': 'sin_impulse'
    },
    chemkin_otomo2018_simple = {
        'mechanism': 'chemkin_otomo2018',
        'R_E': {'type': 'LinearRange', 'start': 1e-06, 'end': 50e-6, 'num_steps': 20},
        'ratio': {'type': 'Const', 'value': 1.0},
        'species': ['H2', 'N2'],
        'fractions': [0.75, 0.25],
        'P_amb': {'type': 'Const', 'value': 1e5},
        'T_inf': {'type': 'Const', 'value': 293.15},
        'alfa_M': {'type': 'Const', 'value': 0.35},
        'P_v': {'type': 'Const', 'value': 2338.1},
        'mu_L': {'type': 'Const', 'value': 0.001},
        'rho_L': {'type': 'Const', 'value': 998.2},
        'c_L': {'type': 'Const', 'value': 1483.0},
        'surfactant': {'type': 'Const', 'value': 1.0},
        'enable_heat_transfer': True,
        'enable_evaporation': True,
        'enable_reactions': True,
        'enable_dissipated_energy': True,
        'target_specie': 'NH3',
        'excitation_params': [
            {'type': 'LinearRange', 'start': -12e5, 'end': -20e5, 'num_steps': 5},
            {'type': 'LinearRange', 'start': 10000.0, 'end': 30000.0, 'num_steps': 5},
            {'type': 'Const', 'value': 1.0}
        ],
        'excitation_type': 'sin_impulse'
    },
    chemkin_kaust2023_n2_simple = {
        'mechanism': 'chemkin_kaust2023_n2',
        'R_E': {'type': 'LinearRange', 'start': 1e-06, 'end': 100e-6, 'num_steps': 20},
        'ratio': {'type': 'Const', 'value': 1.0},
        'species': ['H2', 'N2'],
        'fractions': [0.75, 0.25],
        'P_amb': {'type': 'Const', 'value': 1e5},
        'T_inf': {'type': 'Const', 'value': 293.15},
        'alfa_M': {'type': 'Const', 'value': 0.35},
        'P_v': {'type': 'Const', 'value': 2338.1},
        'mu_L': {'type': 'Const', 'value': 0.001},
        'rho_L': {'type': 'Const', 'value': 998.2},
        'c_L': {'type': 'Const', 'value': 1483.0},
        'surfactant': {'type': 'Const', 'value': 1.0},
        'enable_heat_transfer': True,
        'enable_evaporation': True,
        'enable_reactions': True,
        'enable_dissipated_energy': True,
        'target_specie': 'NH3',
        'excitation_params': [
            {'type': 'LinearRange', 'start': -15e5, 'end': -25e5, 'num_steps': 5},
            {'type': 'LinearRange', 'start': 10000.0, 'end': 80000.0, 'num_steps': 5},
            {'type': 'Const', 'value': 3.0}
        ],
        'excitation_type': 'sin_impulse'
    },
    chemkin_otomo2018_without_O2 = {
        'mechanism': 'chemkin_otomo2018_without_o',
        'R_E': {'type': 'LinearRange', 'start': 0.5e-06, 'end': 10e-6, 'num_steps': 20},
        'ratio': {'type': 'Const', 'value': 1.0},
        'species': ['H2', 'N2'],
        'fractions': [0.75, 0.25],
        'P_amb': {'type': 'Const', 'value': 0.1e5},
        'T_inf': {'type': 'Const', 'value': 293.15},
        'alfa_M': {'type': 'Const', 'value': 0.35},
        'P_v': {'type': 'Const', 'value': 2338.1},
        'mu_L': {'type': 'Const', 'value': 0.001},
        'rho_L': {'type': 'Const', 'value': 998.2},
        'c_L': {'type': 'Const', 'value': 483.0},
        'surfactant': {'type': 'Const', 'value': 0.05},
        'enable_heat_transfer': True,
        'enable_evaporation': False,
        'enable_reactions': True,
        'enable_dissipated_energy': True,
        'target_specie': 'NH3',
        'excitation_params': [
            {'type': 'LinearRange', 'start': -20e5, 'end': -30e5, 'num_steps': 5},
            {'type': 'LinearRange', 'start': 5000.0, 'end': 10000.0, 'num_steps': 5},
            {'type': 'Const', 'value': 5.0}
        ],
        'excitation_type': 'sin_impulse'
    },
    chemkin_kaust2023_n2_high_T = {
        'mechanism': 'chemkin_kaust2023_n2',
        'R_E': {'type': 'LinearRange', 'start': 10e-06, 'end': 100e-6, 'num_steps': 15},
        'ratio': {'type': 'Const', 'value': 1.0},
        'species': ['H2', 'O2', 'N2', 'AR', 'HE'],
        'fractions': [0.45, 0.1, 0.15, 0.2, 0.1],
        'P_amb': {'type': 'Const', 'value': 5e5},
        'T_inf':  {'type': 'LinearRange', 'start': 693.15, 'end': 2093.15, 'num_steps': 10},
        'alfa_M': {'type': 'Const', 'value': 0.35},
        'P_v': {'type': 'Const', 'value': 2338.1},
        'mu_L': {'type': 'Const', 'value': 0.001},
        'rho_L': {'type': 'Const', 'value': 998.2},
        'c_L': {'type': 'Const', 'value': 1483.0},
        'surfactant': {'type': 'Const', 'value': 1.0},
        'enable_heat_transfer': False,
        'enable_evaporation': False,
        'enable_reactions': True,
        'enable_dissipated_energy': False,
        'target_specie': 'NH3',
        'excitation_params': [
            {'type': 'Const', 'value': -20e5},
            {'type': 'Const', 'value': 25000.0},
            {'type': 'Const', 'value': 1.0}
        ],
        'excitation_type': 'sin_impulse'
    },
)

def save_results_to_dict(all_data: pd.DataFrame, study_dir: str, runtime: float) -> dict:
    good_data = all_data.loc[all_data['success'] == True]
    best_energy_demand = good_data['energy_demand'].min()
    avg_runtime = good_data['runtime'].mean()
    avg_num_steps = good_data['num_steps'].mean()
    num_timeout = len(all_data.loc[all_data['type'] == 'timeout'])
    num_semisuccess = len(all_data.loc[all_data['t_last'] > semi_success_time])
    num_success = len(good_data)
    num_total = len(all_data)

    return dict(
        total_runtime = runtime,
        study_dir = study_dir,
        best_energy_demand = best_energy_demand,
        avg_runtime = avg_runtime,
        avg_num_steps = avg_num_steps,
        num_timeout = num_timeout,
        num_semisuccess = num_semisuccess,
        num_success = num_success,
        num_total = num_total,
        rate = 100*num_success / num_total,
    )

def main():
    # Parse command-line arguments
    if len(sys.argv) != 3:
        print("Usage: python3 ./dev/stability_test.py <txtname> <detailed_description>")
        sys.exit(1)

    txtname = sys.argv[1]
    detailed_description = sys.argv[2]

    # Ensure base_dir exists and purge old files except .txt files
    if os.path.exists(base_dir):
        for item in os.listdir(base_dir):
            item_path = os.path.join(base_dir, item)
            if os.path.isdir(item_path):
                shutil.rmtree(item_path)
    else:
        os.makedirs(base_dir)

    # Create summary txt file
    output_file_path = os.path.join(base_dir, f"{txtname}.txt")
    with open(output_file_path, 'a') as output_file:
        output_file.write(f"{txtname}.txt\n")
        output_file.write(f"{detailed_description}\n\n")

    results = dict()
    combined_data = pd.DataFrame()
    start = time.time()

    # Run parameter studies
    for name, parameter_study in parameter_studies.items():
        study_dir = os.path.join(base_dir, name)

        # Run current parameter study
        start_loc = time.time()
        api.run_parameter_study(
            parameter_study=parameter_study,
            t_max = 1.0,
            timeout = timeout,
            save_directory=study_dir
        )
        end_loc = time.time()
        study_dir = study_dir + '1'

        # Load results
        all_data = api.read_parameter_study(study_dir)
        combined_data = pd.concat([combined_data, all_data], ignore_index=True)
        results[name] = save_results_to_dict(all_data, study_dir, end_loc-start_loc)
        with open(output_file_path, 'a') as output_file:
            output_file.write(f"{results[name]}\n")
        
    # Combine results
    end = time.time()
    results['combined'] = save_results_to_dict(combined_data, base_dir, end-start)
    
    # Convert results to a pandas DataFrame
    print('\n')
    results_df = pd.DataFrame.from_dict(results, orient='index')
    results_df = results_df[['total_runtime', 'best_energy_demand', 'avg_runtime', 'avg_num_steps', 'num_timeout', 'num_semisuccess', 'num_success', 'num_total', 'rate']]
    results_df['total_runtime'] = results_df['total_runtime'].map('{:.1f} [s]'.format)
    results_df['best_energy_demand'] = results_df['best_energy_demand'].map('{:.1f} [MJ/kg]'.format)
    results_df['avg_runtime'] = results_df['avg_runtime'].map('{:.3f} [s]'.format)
    results_df['avg_num_steps'] = results_df['avg_num_steps'].map('{:.0f}'.format)
    results_df['rate'] = results_df['rate'].map('{:.2f} [%]'.format)

    print(results_df.to_string(index=True, justify='center'))
    with open(output_file_path, 'a') as output_file:
        output_file.write('\n\n')
        output_file.write(results_df.to_string(index=True, justify='center'))

if __name__ == "__main__":
    main()
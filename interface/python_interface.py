"""
A simple JSON based interface to use Bubble_dynamics_simulation_CPP from Python.
Usage: 
```
from interface import python_interface as api
help(api)
```

You may run the C++ simulation in two modes:
 1. Run a single simulation with the provided control parameters:
    ```
    cpar = api.example_cpar()
    data = api.run_simulation(cpar)
    ```
    You may modify the cpar dictionary or create your own. Results are returned in the data dictionary.
    This include the original control parameters ('cpar') and the numerical solution ('sol'), along with some other info. 
 2. Run a parameter study with the provided parameter study dictionary. (run_parameter_study())
    ```
    parameter_study = api.example_parameter_study()
    api.run_parameter_study(parameter_study, save_directory='./_parameter_studies/test')
    ```
    You may modify the parameter_study dictionary or create your own. Results are saved in a directory.

Afterward you may read results of a parameter study with the read_parameter_study() function:
```
all_data = api.read_parameter_study('./_parameter_studies/test')
good_data = all_data.loc[all_data['success'] == True]
```
You may also convert a line of the all_data DataFrame to a dictionary with the line_to_dict() function:
```
cpar = api.line_to_dict(all_data.iloc[0])
data = api.run_simulation(cpar)
api.plot(data)
```
"""

import numpy as np
import pandas as pd
from scipy.signal import argrelmin
import matplotlib.pyplot as plt
import json
import os
from subprocess import Popen, PIPE, STDOUT
from threading import Thread
import time


def example_cpar() -> dict:
    """
    Returns an example control parameter dictionary.
    Use: cpar = example_cpar()
    """

    file_name = 'interface/example_cpar.json'
    with open(file_name, 'r') as f:
        cpar = json.load(f).get('cpar', {})
    return cpar


def _read_json_and_binary(file_path: str) -> dict:
    """
    Reads Bubble_dynamics_simulation_CPP's output file, which is a JSON binary mix. It begins with an utf8 based JSON string
    and then contains binary data after <BINARY>. The JSON part is parsed, and the binary data is read into numpy arrays.

    Args:
        file_path (str): Path to the file.

    Returns:
        dict: A dictionary containing the JSON data and the binary arrays.
    """

    with open(file_path, 'rb') as f:
        # Read the file until the `<BINARY>` marker
        content = f.read()
        binary_marker = b"<BINARY>"
        marker_index = content.find(binary_marker)
        
        if marker_index == -1:
            raise ValueError("The file does not contain the '<BINARY>' marker.")

        # Split the JSON and binary parts
        json_part = content[:marker_index].decode('utf-8').strip()
        binary_part = content[marker_index + len(binary_marker):]

        # Parse the JSON part
        data = json.loads(json_part)

        # Extract metadata from the JSON
        if "sol" not in data:
            raise KeyError("The JSON part does not contain the 'sol' key.")
        if "num_saved_steps" not in data["sol"]:
            raise KeyError("The 'sol' key does not contain the 'num_saved_steps' key.")
        if "num_dim" not in data["sol"]:
            raise KeyError("The 'sol' key does not contain the 'num_dim' key.")
        num_saved_steps = data["sol"]["num_saved_steps"]
        num_dim = data["sol"]["num_dim"]

        # Read the binary part
        t = np.frombuffer(binary_part[:num_saved_steps * 8], dtype=np.float64)
        x = np.frombuffer(binary_part[num_saved_steps * 8:], dtype=np.float64).reshape((num_saved_steps, num_dim))

        # Add the binary data to the dictionary
        data = dict(data)
        data["sol"]["t"] = t
        data["sol"]["x"] = x
        data["sol"]["total_error"] = np.array(data["sol"]["total_error"])

    return data


def _check_path(
        path: str
) -> str:
    """
    Checks if the path exists and returns the absolute path.
    """

    if not os.path.exists(path):
        raise FileNotFoundError(f'The path "{path}" does not exist.')
    else:
        path = os.path.abspath(path)
    return path


def _run_cpp_simulation(
    command_list: list[str]
) -> int:
    """
    Runs the C++ simulation with the given command list. Errors and outputs are printed to the console.
    Returns the return code of the process.
    """


    # Functions to read stdout and stderr
    def stream_reader(stream, prefix):
        for line in iter(stream.readline, ''):
            print(f'{prefix} {line}', end='')
        stream.close()

    with Popen(command_list, stdout=PIPE, stderr=PIPE, text=True) as p:
        # Create threads to read stdout and stderr concurrently
        stdout_thread = Thread(target=stream_reader, args=(p.stdout, '[stdout]'))
        stderr_thread = Thread(target=stream_reader, args=(p.stderr, '[stderr]'))

        # Start the threads
        stdout_thread.start()
        stderr_thread.start()

        # Wait for the threads to finish
        stdout_thread.join()
        stderr_thread.join()

        # Wait for the process to complete and get the return code
        return_code = p.wait()
        return return_code
    

def run_simulation(
    cpar: dict,
    json_path: str = 'ignore.json',
    executable_path: str = './bin/main' if os.name == 'posix' else './bin/main.exe',
    t_max: float = 1.0,
    timeout: float = 60.0,
    save_steps: bool = True,
) -> dict:
    """
    Runs the Bubble_dynamics_simulation_CPP executable in --run mode to run a single simulation.
    Saves the provided control parameters dictionary to a JSON file and uses it as input to run the simulation.
    Captures stdout and stderr in real-time and prints them as they arrive.

    Args:
        cpar (dict): Parameters for the simulation.
        json_path (str): Path to the JSON file.
        executable_path (str): Path to the executable.
        t_max (float): Maximum simulation time.
        timeout (float): Timeout for the simulation.
        save_steps (bool): Whether to save the simulation steps.

    Returns:
        dict: A dictionary containing the simulation results.
              See key 'sol' for the numerical solution and 'cpar' for the original control parameters.
    """
    
    # Write the JSON file
    json_str = json.dumps(dict(cpar=cpar), indent=4)
    with open(json_path, 'w') as f:
        f.write(json_str)
    
    # Check paths
    json_path = _check_path(json_path)
    executable_path = _check_path(executable_path)

    # Run simulation (call the executable)
    start = time.time()
    command_list = [
        executable_path, '--run', json_path,
        '--tmax', str(t_max),
        '--timeout', str(timeout),
    ]
    if save_steps: command_list.append('--save')

    return_code = _run_cpp_simulation(command_list)
    print(f'\n{executable_path} returned with code {return_code} after {time.time() - start:.4f} seconds.')

    # Read the JSON-binary file
    data = _read_json_and_binary(json_path)
    return data


def example_parameter_study() -> dict:
    """
    Returns an example dictionary for a parameter study.
    Usage: parameter_study = example_parameter_study()
    Types:
     * {"type": "Const", value=1.0}
     * {"type": "LinearRange", start=0.0, end=1.0, num_steps=10}
     * {"type": "PowRange", start=1.0, end=100.0, num_steps=10, base=2.0}
    """

    file_name = 'interface/example_parameter_study.json'
    with open(file_name, 'r') as f:
        parameter_study = json.load(f).get('parameter_study', {})
    return parameter_study


def run_parameter_study(
    parameter_study: dict,
    json_path: str = 'ignore.json',
    executable_path: str = './bin/main' if os.name == 'posix' else './bin/main.exe',
    t_max: float = 1.0,
    timeout: float = 60.0,
    save_directory: str = './_parameter_studies/test',
):
    """
    Runs the Bubble_dynamics_simulation_CPP executable in --parameter_study mode to run a bruteforce parameter study.
    Saves the provided parameter_study dictionary to a JSON file and uses it as input to run the study.
    Captures stdout and stderr in real-time and prints them as they arrive.

    Args:
        parameter_study (dict): Parameters for the simulation.
        json_path (str): Path to the JSON file.
        executable_path (str): Path to the executable.
        t_max (float): Maximum simulation time.
        timeout (float): Timeout for the simulation.
        save_directory (str): Where to save the parameter study.
    """
    
    # Write the JSON file
    json_str = json.dumps(dict(parameter_study=parameter_study), indent=4)
    with open(json_path, 'w') as f:
        f.write(json_str)
    
    # Check paths
    json_path = _check_path(json_path)
    executable_path = _check_path(executable_path)

    # Run simulation (call the executable)
    start = time.time()
    command_list = [
        executable_path, '--parameter_study', json_path,
        '--tmax', str(t_max),
        '--timeout', str(timeout),
        '--directory', str(save_directory),
    ]

    return_code = _run_cpp_simulation(command_list)
    print(f'\n{executable_path} returned with code {return_code} after {time.time() - start:.4f} seconds.')


def _print_data(data, print_it=True):
    """
    Prints the data dictionary in an organized way.
    Arguments:
     * data: data dictionary to be printed
     * print_it: if True, the function will print the text. If False, it will return the text (string)
    """

    cpar = data.get('cpar', {})
    sol = data.get('sol', {})
    
    # Control parameters
    text = "Control Parameters:\n"
    text += f"  ID: {cpar.get('ID', 'N/A')}\n"
    text += f"  Mechanism: {cpar.get('mechanism', 'N/A')}\n"
    text += f"  R_E: {1e6 * cpar.get('R_E', 'N/A')} [um]\n"
    text += f"  Species: {cpar.get('species', [])}\n"
    text += f"  Fractions: {cpar.get('fractions', [])}\n"
    text += f"  P_amb: {cpar.get('P_amb', 'N/A')} [Pa]\n"
    text += f"  T_inf: {cpar.get('T_inf', 'N/A')} [K]\n"
    text += f"  alfa_M: {cpar.get('alfa_M', 'N/A')} [-]\n"
    text += f"  P_v: {cpar.get('P_v', 'N/A')} [Pa]\n"
    text += f"  mu_L: {cpar.get('mu_L', 'N/A')} [Pa·s]\n"
    text += f"  rho_L: {cpar.get('rho_L', 'N/A')} [kg/m³]\n"
    text += f"  c_L: {cpar.get('c_L', 'N/A')} [m/s]\n"
    text += f"  Surfactant: {cpar.get('surfactant', 'N/A')}\n"
    text += f"  Enable Heat Transfer: {cpar.get('enable_heat_transfer', 'N/A')}\n"
    text += f"  Enable Evaporation: {cpar.get('enable_evaporation', 'N/A')}\n"
    text += f"  Enable Reactions: {cpar.get('enable_reactions', 'N/A')}\n"
    text += f"  Enable Dissipated Energy: {cpar.get('enable_dissipated_energy', 'N/A')}\n"
    text += f"  Target Specie: {cpar.get('target_specie', 'N/A')}\n"
    text += f"  Excitation Params: {cpar.get('excitation_params', 'N/A')}\n"
    text += f"  Excitation Type: {cpar.get('excitation_type', 'N/A')}\n"

    # Simulation info
    text += "\nSimulation Info:\n"
    text += f"  Success: {sol.get('success', 'N/A')}\n"
    text += f"  Error: {sol.get('error', 'N/A')}\n"
    text += f"  Runtime: {sol.get('runtime', 'N/A')} [s]\n"
    text += f"  Num Steps: {sol.get('num_steps', 'N/A')}\n"
    text += f"  Num Repeats: {sol.get('num_repeats', 'N/A')}\n"
    text += f"  Num Function Evaluations: {sol.get('num_fun_evals', 'N/A')}\n"
    text += f"  Num Jacobian Evaluations: {sol.get('num_jac_evals', 'N/A')}\n"
    text += f"  t_last = {sol.get('t', ['N/A'])[-1]} [s]\n"

    # Results
    text += "\nResults:\n"
    text += f"  Dissipated Energy: {data.get('dissipated_energy', 'N/A')} [J]\n"
    text += f"  n_target_specie: {data.get('n_target_specie', 'N/A')} [mol]\n"
    text += f"  Energy Demand: {data.get('energy_demand', 'N/A')} [MJ/kg]\n"

    if print_it:
        print(text)
    else:
        return text
    

def plot(data, n=5.0, base_name='', format='png',
         presentation_mode=False, show_legend=False, show_cpar=True):
    """
    This funfction plots the results of the simulation form data (returned by run_simulation()).
    Parameters:
     * data: data dictionary (returned by run_simulation()) | 
     * n: how long should the plotted time interval be compared to the collapse time (default: 5 [-])
     * base_name: save plots as image (default: '' alias do not save) | 
               use base_name='plot' --> plot_1.png, plot_2.png |  
               use base_name='images/plot' to save into images folder |  
               using a folder for images is recommend |  
               this folder have to be created manually
     * format: format of the saved images (available: png, pdf, ps, eps, svg)
     * presentation_mode: if True, the plot will be in presentation mode (default: False)
     * show_legend: if True, the legend will be visible with every single species (default: False)
     * show_cpar: if True, the control parameters will be printed on the plot (default: False)
    """

# Calculations 
    cpar = data['cpar']
    sol = data['sol']
    x = sol['x']
    t = sol['t']
    loc_min = argrelmin(x[:, 0])
    if not len(loc_min) == 0 and not len(loc_min[0]) == 0:
        collapse_time = t[loc_min[0][0]] # collapse time (first loc min of R) [s]
    else:
        collapse_time = t[-1]

    t_last = n * collapse_time
    if t_last < 1e-7 or t[-1] < t_last or n < 0 or not sol['success']:
        end_index = -1
    else:
        end_index = np.where(t >= t_last)[0][0]

    if t[end_index] < 1e-3:
        t = t[:end_index] * 1e6 # [us]
        t_label = '$t$ [μs]'
    else:
        t = t[:end_index] * 1e3 # [ms]
        t_label = '$t$ [ms]'
    R = x[:end_index, 0] # [m]
    R_dot = x[:end_index, 1] # [m/s]
    T = x[:end_index, 2] # [K]
    c = x[:end_index, 3:-1] # [mol/cm^3]

    V = 4.0 / 3.0 * (100.0 * R) ** 3 * np.pi # [cm^3]
    n = (c.T * V).T
    #internal_pressure = np.sum(n, axis=1) * 8.31446 * T / V # [MPa]

# plot R and T
    linewidth = 2.0 if presentation_mode else 1.0
    plt.rcParams.update({'font.size': 24 if presentation_mode else 18})
    fig1 = plt.figure(figsize=(16, 9) if presentation_mode else (20, 6))
    ax1 = fig1.add_subplot(axisbelow=True)
    ax2 = ax1.twinx() 
    ax1.plot(t, R / cpar['R_E'], color = 'b', linewidth = linewidth)
    ax2.plot(t, T, color = 'r', linewidth = linewidth, linestyle = '-.')

    ax1.set_xlabel(t_label)
    ax1.set_ylabel('$R/R_E$ [-]', color = 'b')
    ax2.set_ylabel('$T$ [K]', color = 'r')
    if not presentation_mode: ax1.grid()
    
# textbox with initial conditions
    text = f'Initial conditions:\n'
    text += f'    $R_E$ = {1e6*cpar["R_E"]: .2f} $[\mu m]$\n'
    #if cpar['ratio'] != 1.0:
    #    text += f'    $R_0/R_E$ = {cpar['ratio']: .2f} $[-]$\n'
    text += f'    $P_{{amb}}$ = {1e-5*cpar["P_amb"]: .2f} $[bar]$\n'
    text += f'    $T_{{inf}}$ = {cpar["T_inf"]-273.15: .2f} $[°C]$\n'
    text += f'    $P_{{vapour}}$ = {cpar["P_v"]: .1f} $[Pa]$\n'
    text += f'Initial content:\n    '
    for gas, fraction in zip(cpar['species'], cpar['fractions']):
        text += f'{int(100*fraction)}% {gas}, ' 
    text = text[:-2] + f'\nExcitation = {cpar["excitation_type"]}:\n'
    params = cpar['excitation_params']
    names = data['excitation']['names']
    units = data['excitation']['units']
    for i, value in enumerate(params):
        text += f'    {names[i]} = {value} [{units[i]}],\n'
    text = text[:-2]

    if show_cpar and not presentation_mode:
        ax2.text(
            0.98, 0.95, # coordinates
            text, transform=ax1.transAxes,
            horizontalalignment='right', verticalalignment='top', multialignment='left',
            fontsize=14, fontstyle='oblique',
            bbox={'facecolor': 'white', 'alpha': 1.0, 'pad': 10},
        )
    
    plt.show()

# plot reactions
    plt.rcParams.update({'font.size': 24 if presentation_mode else 18})
    fig2 = plt.figure(figsize=(16, 9) if presentation_mode else (20, 9))
    ax = fig2.add_subplot(axisbelow=True)

    # plot the lines
        # use this to generate colors:
            # import seaborn as sns
            # colors = sns.color_palette('Set1', n_colors=10)
            # print(colors.as_hex()); colors
    colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#d48044', '#33adff', '#a65628', '#f781bf', '#d444ca', '#d4ae44']
    color_index = 0
    texts = []
    max_mol = np.max(np.nan_to_num(n, nan=0.0, posinf=0.0, neginf=0.0), axis=1) # maximum amounts of species [mol]
    indexes_to_plot = []
    for i, specie in enumerate(data['mechanism']['species_names']):
        if n[-1, i] > 1e-24:
            indexes_to_plot.append(i)
    for i, specie in enumerate(data['mechanism']['species_names']):
        name = specie
        for digit in range(10): # turns 'H2O2' into 'H_2O_2'
            name = name.replace(str(digit), '_' + str(digit))
        if i in indexes_to_plot:
            color = colors[color_index]
            color_index = color_index + 1 if color_index < len(colors) - 1 else 0
            linewidth = 2.0 if presentation_mode and n[-1, i] > 1e-24 else 1.0
            ax.plot(t, n[:, i], linewidth = linewidth, color=color, label = '$' + name + '$') # PLOT HERE
            texts.append((color, name, n[-1, i]))            
        elif not presentation_mode:
            linewidth = 2.0 if presentation_mode else 1.0
            ax.plot(t, n[:, i], linewidth = linewidth, label = '$' + name + '$')  # PLOT HERE

    # make legend
    texts.sort(key=lambda x: x[-1], reverse=True)
    last_n_last = 1.0e100
    for text in texts:
        color, name, n_last = text
        # spaceing
        if n_last < 1e-24: continue
        limit = 5.5 if presentation_mode else 3.5
        if last_n_last / n_last < limit:
            n_last = last_n_last / limit
        last_n_last = n_last
        # place text
        ax.text(
            t[-1],
            n_last,
            '$' + name + '$',
            color=color,
            fontsize=24 if presentation_mode else 18,
            verticalalignment='center',
            bbox={'facecolor': 'white', 'pad': 0, 'linewidth': 0.0},
        )

    # plot settings
    ax.set_ylim([1e-24, 5.0*max(max_mol)])
    ax.set_yscale('log')
    ax.set_xlabel(t_label)
    ax.set_ylabel('$n_k$ [mol]')
    if not presentation_mode: ax.grid()
    if show_legend: ax.legend()

    plt.show()
    
# saving the plots
    if base_name != '':
        if format not in ['png', 'pdf', 'ps', 'eps', 'svg']:
            print(f'Invalid image format {format}, png is used instead. ')
            format = 'png'
        name = base_name + '_radial.' + format
        if os.path.exists(name):
            print(f'Error: {name} already exists. ')
            return None
        try:
            if format == 'png':
                metadata = {key: str(data[key]) for key in data.keys()}
            else:
                metadata = {}
            fig1.savefig(base_name+'_radial.'+format, format=format, metadata=metadata, bbox_inches='tight')
            fig2.savefig(base_name+'_molar.'+format, format=format, metadata=metadata, bbox_inches='tight')
        except Exception as error:
            print(f'Error in saving {base_name}_1.png')
            print(error)

# print data
    _print_data(data, print_it=True)
    return None


def read_parameter_study(directory: str) -> pd.DataFrame:
    """
    Read a parameter study directory and return a pandas DataFrame.
    Usage:
        all_data = read_parameter_study('path/to/parameter/study')
        good_data = all_data.loc[all_data['success'] == True]

    You may also retrive the settings:
        path = 'path/to/parameter/study/bruteforce_parameter_study_settings.json'
        with open(path, 'r') as f:
            settings = json.load(f)
    """

    # create a dataframe
    print(f'Found files:')
    directory = _check_path(directory)
    all_data = pd.DataFrame()
    num = 0

    # iterate trough all files in directory (including subdirectories)
    for (root, dirs, files) in os.walk(directory):
        for file in files:
            if 'ipynb_checkpoints' in root: # ignore python rubish
                continue
            if file[-4:] != '.csv':
                continue

            # read file
            num += 1
            current_data = pd.read_csv(os.path.join(root, file))

            # Cast object-dtype columns with all-bool values to bool dtype
            for col in current_data.columns:
                if current_data[col].dtype == 'object' and all(current_data[col].dropna().map(lambda x: isinstance(x, bool))):
                    current_data[col] = current_data[col].astype(bool)

            subdir = os.path.join(root.removeprefix(directory), file)
            print(f'\t{subdir: <64} ({current_data.shape[0]: >4} rows)')
            all_data = pd.concat([all_data, current_data])
        
    # Print some stats:
    print(f'_______________________________________')
    total = all_data.shape[0]
    all_data = all_data.sort_values(['energy_demand'], ascending=True)
    print(f'total: {total: >4} rows in {num: >2} files')

    return all_data


def line_to_dict(line):
    """
    Convert a line of all_data (from read_parameter_study()) to a dictionary containing the control parameters (cpar).
    Usage: cpar = line_to_dict(all_data.iloc[0])
    """

    species = str(line['species'])
    species = [str(specie) for specie in species.split(';') if specie != '']
    fractions = str(line['fractions'])
    fractions = [float(frac) for frac in fractions.split(';') if frac != '']
    excitation_params = [float(param) for param in str(line['excitation_params']).split(';') if param != '']

    return dict(
        ID = int(line['ID']),
        mechanism = str(line['mechanism']),
        R_E = float(line['R_E']),
        #ratio = line['ratio'],
        species = species,
        fractions = fractions,
        P_amb = float(line['P_amb']),
        T_inf = float(line['T_inf']),
        alfa_M = float(line['alfa_M']),
        P_v = float(line['P_v']),
        mu_L = float(line['mu_L']),
        rho_L = float(line['rho_L']),
        c_L = float(line['c_L']),
        surfactant = float(line['surfactant']),
        enable_heat_transfer = bool(line['enable_heat_transfer']),
        enable_evaporation = bool(line['enable_evaporation']),
        enable_reactions = bool(line['enable_reactions']),
        enable_dissipated_energy = bool(line['enable_dissipated_energy']),
        target_specie = str(line['target_specie']),
        excitation_params = excitation_params,
        excitation_type = str(line['excitation_type']),
    )
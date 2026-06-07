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
import matplotlib.pyplot as plt
import json
import os
import shutil
from subprocess import Popen, PIPE, STDOUT
from threading import Thread
import time

# Loading the interface
system = None
executable_path = None
if os.name == 'nt':
    if os.path.exists('./bin/main.exe'):
        system = 'windows'
        executable_path = './bin/main.exe'
    elif os.path.exists('./bin/main'):
        system = 'wsl'
        executable_path = './bin/main'
elif os.name == 'posix':
    if os.path.exists('./bin/main'):
        system = 'linux'
        executable_path = './bin/main'

print(f'{executable_path=}, {system=}')
if executable_path is None:
    print('Warning: No default executable path found. Please check the ./bin directory. Perhaps you forgot to compile with ./dev/build.py')
if system is None:
    print('Warning: System is not supported. You may use:')
    print('    Windows with ./bin/main.exe')
    print('    WSL with ./bin/main')
    print('    Linux with ./bin/main')


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
    
    The output dictionary includes:
        - t: time array
        - R, R_dot, T, E_diss: individual state variables
        - p_excitation, p_internal: pressure time series
        - n_<species_name>: individual species molar amounts [mol]
        - x: full state array (for backward compatibility)

    Args:
        file_path (str): Path to the file.

    Returns:
        dict: A dictionary containing the JSON data and the binary arrays with restructured sol data.
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

        # Read the binary part sequentially:
        offset = 0
        t = np.frombuffer(binary_part[offset:offset + num_saved_steps * 8], dtype=np.float64)
        offset += num_saved_steps * 8
        
        x = np.frombuffer(binary_part[offset:offset + num_saved_steps * num_dim * 8], dtype=np.float64).reshape((num_saved_steps, num_dim))
        offset += num_saved_steps * num_dim * 8
        
        p_excitation = np.frombuffer(binary_part[offset:offset + num_saved_steps * 8], dtype=np.float64)
        offset += num_saved_steps * 8
        
        p_internal = np.frombuffer(binary_part[offset:offset + num_saved_steps * 8], dtype=np.float64)

        # Add the binary data to the dictionary and restructure sol
        data = dict(data)
        sol = dict(data["sol"])
        
        sol["t"] = t
        sol["x"] = x  # Keep for backward compatibility
        sol["p_excitation"] = p_excitation
        sol["p_internal"] = p_internal
        sol["total_error"] = np.array(sol.get("total_error", []))
        
        sol["R"] = x[:, 0]           # radius [m]
        sol["R_dot"] = x[:, 1]       # radius velocity [m/s]
        sol["T"] = x[:, 2]           # temperature [K]
        sol["E_diss"] = x[:, -1]     # dissipated energy [J]
        
        if "mechanism" in data and "species_names" in data["mechanism"]:
            species_names = data["mechanism"]["species_names"]
            concentrations = x[:, 3:-1]
            V = (4.0 / 3.0) * np.pi * x[:, 0]**3  # bubble volume [m^3]
            for i, species_name in enumerate(species_names):
                sol[f"n_{species_name}"] = concentrations[:, i] * V  # molar amount [mol]
        
        data["sol"] = sol

    return data


def _check_path(
        path: str
) -> str:
    """
    Checks if the path exists and returns the absolute path.
    Does nothing if system is 'wsl' (Windows Subsystem for Linux).
    """

    if system == 'wsl':
        return path # no checks

    if not os.path.exists(path):
        raise FileNotFoundError(f'The path "{path}" does not exist.')
    else:
        path = os.path.abspath(path).replace('\\\\', '/').replace('\\', '/')
    return path


def _run_cpp_simulation(
    command_list: list[str],
    show_stdout: bool = True,
    show_stderr: bool = True
) -> int:
    """
    Runs the C++ simulation with the given command list. Errors and outputs are printed to the console.
    Returns the return code of the process.
    """

    if system == 'wsl':
        command_list = ['wsl'] + command_list

    # Functions to read stdout and stderr
    def stream_reader(stream, prefix, show_output):
        for line in iter(stream.readline, ''):
            if show_output:
                print(f'{prefix} {line}', end='')
        stream.close()

    with Popen(command_list, stdout=PIPE, stderr=PIPE, text=True) as p:
        # Create threads to read stdout and stderr concurrently
        stdout_thread = Thread(target=stream_reader, args=(p.stdout, '[stdout]', show_stdout))
        stderr_thread = Thread(target=stream_reader, args=(p.stderr, '[stderr]', show_stderr))

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
    executable_path: str = executable_path,
    t_max: float = 1.0,
    timeout: float = 60.0,
    save_steps: bool = True,
    save_jacobian: bool = False,
    show_output: bool = True
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
        save_jacobian (bool): Whether to save the Jacobian matrixes.
        show_output (bool): Whether to show the stdout and stderr of the simulation.

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
    if save_jacobian and os.path.exists('./_jacobians'):
        shutil.rmtree('./_jacobians')

    # Run simulation (call the executable)
    start = time.time()
    command_list = [executable_path,
        '--run', json_path,
        '--tmax', str(t_max),
        '--timeout', str(timeout),
    ]
    if save_steps: command_list.append('--save')
    if save_jacobian: command_list.append('--save_jacobian')

    return_code = _run_cpp_simulation(command_list, show_stdout=show_output, show_stderr=show_output)
    if show_output:
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
    executable_path: str = executable_path,
    t_max: float = 1.0,
    timeout: float = 60.0,
    save_directory: str = './_parameter_studies/test',
    cpu_count: int = None,
    show_stdout: bool = True,
    show_stderr: bool = True
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
        cpu_count (int): Number of CPUs to use for the simulation.
        show_stdout (bool): Whether to show the stdout of the simulation.
        show_stderr (bool): Whether to show the stderr of the simulation.
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
    if cpu_count is not None:
        command_list.append('--cpu')
        command_list.append(str(cpu_count))

    return_code = _run_cpp_simulation(command_list, show_stdout=show_stdout, show_stderr=show_stderr)
    print(f'\n{executable_path} returned with code {return_code} after {time.time() - start:.4f} seconds.')


def _print_data(data, print_it=True):
    """
    Prints the data dictionary in an organized way.
    Arguments:
     * data: data dictionary to be printed
     * print_it: if True, the function will print the text. If False, it will return the text (string)
    """

    nan = float('nan')
    cpar = data.get('cpar', {})
    sol = data.get('sol', {})
    postproc = data.get('postproc', {})

    # Control parameters
    text = "Control Parameters:\n"
    text += f"  ID: {cpar.get('ID', nan)}\n"
    text += f"  Mechanism: {cpar.get('mechanism', nan)}\n"
    text += f"  R_E: {1e6 * cpar.get('R_E', nan)} [um]\n"
    text += f"  Ratio: {cpar.get('ratio', nan)} [-]\n"
    text += f"  Species: {cpar.get('species', [])}\n"
    text += f"  Fractions: {cpar.get('fractions', [])}\n"
    text += f"  P_amb: {cpar.get('P_amb', nan)} [Pa]\n"
    text += f"  T_inf: {cpar.get('T_inf', nan)} [K]\n"
    text += f"  alpha_M: {cpar.get('alpha_M', nan)} [-]\n"
    text += f"  P_v: {cpar.get('P_v', nan)} [Pa]\n"
    text += f"  mu_L: {cpar.get('mu_L', nan)} [Pa·s]\n"
    text += f"  rho_L: {cpar.get('rho_L', nan)} [kg/m³]\n"
    text += f"  c_L: {cpar.get('c_L', nan)} [m/s]\n"
    text += f"  Surfactant: {cpar.get('surfactant', nan)}\n"
    text += f"  Bubble Dynamics: {cpar.get('bubble_dynamics', nan)}\n"
    text += f"  Liquid EOS Params: {cpar.get('liquid_eos_params', nan)}\n"
    text += f"  Enable Heat Transfer: {cpar.get('enable_heat_transfer', nan)}\n"
    text += f"  Enable Evaporation: {cpar.get('enable_evaporation', nan)}\n"
    text += f"  Enable Reactions: {cpar.get('enable_reactions', nan)}\n"
    text += f"  Enable Dissipated Energy: {cpar.get('enable_dissipated_energy', nan)}\n"
    text += f"  Enable Van der Waals: {cpar.get('enable_van_der_waals', nan)}\n"
    text += f"  Enable Rate Thresholding: {cpar.get('enable_rate_thresholding', nan)}\n"
    text += f"  Target Specie: {cpar.get('target_specie', nan)}\n"
    text += f"  Excitation Type: {cpar.get('excitation_type', nan)}\n"
    text += f"  Excitation Params: {cpar.get('excitation_params', nan)}\n"
    text += f"  Excitation Cycles: {cpar.get('excitation_cycles', nan)} [-]\n"
    text += f"  Ramp Up Cycles: {cpar.get('ramp_up_cycles', nan)} [-]\n"

    # Simulation info
    text += "\nSimulation Info:\n"
    text += f"  Success: {sol.get('success', nan)}\n"
    text += f"  Error: {sol.get('error', nan)}\n"
    text += f"  Runtime: {sol.get('runtime', nan):.3f} [s]\n"
    text += f"  Num Steps: {sol.get('num_steps', nan)}\n"
    text += f"  Num Repeats: {sol.get('num_repeats', nan)}\n"
    text += f"  Num Function Evaluations: {sol.get('num_fun_evals', nan)}\n"
    text += f"  Num Jacobian Evaluations: {sol.get('num_jac_evals', nan)}\n"
    t_arr = sol.get('t', [])
    text += f"  t_last = {(t_arr[-1] if len(t_arr)>0 else nan): .4g} [s]\n"
    R_E   = cpar.get('R_E', nan) or nan
    T_inf_v = cpar.get('T_inf', nan) or nan
    R_max_v = postproc.get('R_max', nan)
    R_min_v = postproc.get('R_min', nan)
    T_max_v = postproc.get('T_max', nan)
    T_min_v = postproc.get('T_min', nan)
    text += f"  R_max = {1e6*R_max_v:.3f} [um]  (R_max/R_E = {R_max_v/R_E:.2f})\n"
    text += f"  R_min = {1e6*R_min_v:.3f} [um]  (R_min/R_E = {R_min_v/R_E:.2f})\n"
    text += f"  T_max = {T_max_v:.2f} [K]  (T_max/T_inf = {T_max_v/T_inf_v:.2f})\n"
    text += f"  T_min = {T_min_v:.2f} [K]  (T_min/T_inf = {T_min_v/T_inf_v:.2f})\n"
    text += f"  t_peak = {1e6*postproc.get('t_peak', nan):.6g} [us]\n"
    text += f"  v_max = {postproc.get('v_max', nan):.6g} [m/s]\n"
    text += f"  p_internal_max = {postproc.get('p_internal_max', nan):.6g} [Pa]\n"
    text += f"  p_internal_min = {postproc.get('p_internal_min', nan):.6g} [Pa]\n"
    text += f"  Ma_max = {postproc.get('Ma_max', nan):.6g} [-]\n"
    if cpar.get('bubble_dynamics', 'keller_miksis') != 'keller_miksis':
        text += f"  T_L_max = {postproc.get('T_L_max', nan):.6g} [K]\n"
        text += f"  c_L_max = {postproc.get('c_L_max', nan):.6g} [m/s]\n"
        text += f"  rho_L_max = {postproc.get('rho_L_max', nan):.6g} [kg/m^3]\n"

    # Results
    text += "\nResults:\n"
    text += f"  Dissipated Energy: {postproc.get('dissipated_energy', nan):.6g} [J]\n"
    text += f"  Expansion Work: {postproc.get('expansion_work', nan):.6g} [J]\n"
    text += f"  n_target_specie: {postproc.get('n_target_specie', nan):.6g} [mol]\n"
    text += f"  Energy Demand: {postproc.get('energy_demand', nan):.6g} [MJ/kg]\n"

    if print_it:
        print(text)
    else:
        return text
    

def plot(data, n=5.0, base_name='', format='png',
         presentation_mode=False, show_legend=False, show_cpar=True, plot_pressure=False):
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
     * plot_pressure: if True, plot external excitation and internal pressure (default: False)
    """

# Calculations 
    cpar = data['cpar']
    sol = data['sol']
    x = sol['x']
    t = sol['t']

    t_last = n * data.get('postproc', {}).get('t_peak', 0.0)
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
    c = x[:end_index, 3:-1] # [mol/m^3]

    V = 4.0 / 3.0 * R**3 * np.pi # [m^3]
    n = (c.T * V).T # [mol]
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
    text += f'    $R_E$ = {1e6*cpar["R_E"]: .2f} $[\\mu m]$\n'
    if cpar['ratio'] != 1.0:
        text += f'    $R_0/R_E$ = {cpar["ratio"]: .2f} $[-]$\n'
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
    
# plot pressure excitation and internal pressure
    if plot_pressure:
        p_excitation = sol['p_excitation'][:end_index] * 1e-3
        p_internal = sol['p_internal'][:end_index] * 1e-6
        
        plt.rcParams.update({'font.size': 24 if presentation_mode else 18})
        linewidth = 2.0 if presentation_mode else 1.0
        fig3 = plt.figure(figsize=(16, 9) if presentation_mode else (20, 6))
        ax1 = fig3.add_subplot(axisbelow=True)
        ax2 = ax1.twinx()
        
        ax1.plot(t, p_excitation, color='darkorange', label='external pressure', linewidth=linewidth)
        ax2.plot(t, p_internal, color='g', label='internal pressure', linewidth=linewidth)

        ax1.set_xlabel(t_label)
        ax1.set_ylabel('Pressure excitation [kPa]', color='darkorange')
        ax2.set_ylabel('Internal pressure [MPa]', color='g')
        ax2.set_ylim([0.5*np.min(p_internal), 2.0*np.max(p_internal)])
        ax2.set_yscale('log')
        if not presentation_mode: ax1.grid()

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

    # treat weird tokens as NaN when reading CSVs
    _na_tokens = ['-nan(ind)', 'nan', 'NaN', 'inf', '-inf', 'Infinity', '-Infinity']

    # iterate trough all files in directory (including subdirectories)
    for (root, dirs, files) in os.walk(directory):
        for file in files:
            if 'ipynb_checkpoints' in root: # ignore python rubish
                continue
            if file[-4:] != '.csv':
                continue

            # read file
            num += 1
            current_data = pd.read_csv(
                os.path.join(root, file),
                na_values=_na_tokens,
                keep_default_na=True
            )

            # Cast object-dtype columns with all-bool values to bool dtype
            for col in current_data.columns:
                if current_data[col].dtype == 'object' and all(current_data[col].dropna().map(lambda x: isinstance(x, bool))):
                    current_data[col] = current_data[col].astype(bool)

            subdir = os.path.join(root.removeprefix(directory), file)
            print(f'\t{subdir: <64} ({current_data.shape[0]: >4} rows)')
            all_data = pd.concat([all_data, current_data], ignore_index=True)
        
    # Print some stats:
    print(f'_______________________________________')
    total = all_data.shape[0]

    # Ensure numeric sort key and keep NaNs last
    if 'energy_demand' in all_data.columns:
        all_data['energy_demand'] = pd.to_numeric(all_data['energy_demand'], errors='coerce')
        all_data = all_data.sort_values(['energy_demand'], ascending=True, na_position='last')

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
    liquid_eos_params = [float(param) for param in str(line.get('liquid_eos_params', '')).split(';') if param != '']

    return dict(
        ID = int(line['ID']),
        mechanism = str(line['mechanism']),
        R_E = float(line['R_E']),
        ratio = line['ratio'],
        species = species,
        fractions = fractions,
        P_amb = float(line['P_amb']),
        T_inf = float(line['T_inf']),
        alpha_M = float(line['alpha_M']),
        P_v = float(line['P_v']),
        mu_L = float(line['mu_L']),
        rho_L = float(line['rho_L']),
        c_L = float(line['c_L']),
        surfactant = float(line['surfactant']),
        bubble_dynamics = str(line.get('bubble_dynamics', 'keller_miksis')),
        liquid_eos_params = liquid_eos_params,
        enable_heat_transfer = bool(line['enable_heat_transfer']),
        enable_evaporation = bool(line['enable_evaporation']),
        enable_reactions = bool(line['enable_reactions']),
        enable_dissipated_energy = bool(line['enable_dissipated_energy']),
        enable_van_der_waals = bool(line['enable_van_der_waals']),
        enable_rate_thresholding = bool(line['enable_rate_thresholding']),
        target_specie = str(line['target_specie']),
        excitation_type = str(line['excitation_type']),
        excitation_params = excitation_params,
        excitation_cycles = int(line['excitation_cycles']),
        ramp_up_cycles = int(line['ramp_up_cycles']),
    )
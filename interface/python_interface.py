"""
A simple JSON based interface to use Bubble_dynamics_simulation_CPP from Python.
Usage: 
```Python
from interface import python_interface as api
help(api)
```
You may run the C++ simulation in two modes:
 1. Run a single simulation with the provided control parameters:
    ```Python
    cpar = api.example_cpar()
    data = api.run_simulation(cpar, 'ignore.json')
    ```
    You may modify the cpar dictionary or create your own. Results are returned in an other dictionary.
    This include the original control parameters ('cpar') and the numerical solution ('sol'), along with some other info. 
 2. Run a parameter study with the provided parameter study dictionary. (run_parameter_study())
    ```Python
    parameter_study = api.example_parameter_study()
    api.run_parameter_study(
        parameter_study,
        json_path='ignore.json'
    )
    ```
    You may modify the parameter_study dictionary or create your own. Results are saved in a directory.

Both functions will create a JSON file with the provided parameters and run the C++ simulation with the appropriate command line arguments.
Further post processing functions are available:


"""

import numpy as np
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

    return dict(
        ID = 0,
        mechanism = 'chemkin_ar_he',
        species = ['O2'],
        fractions = [1.0],
        P_amb = 101325.0,
        T_inf = 293.15,
        alfa_M = 0.35,
        P_v = 2338.1,
        mu_L = 0.001,
        rho_L = 998.2,
        c_L = 1483.0,
        surfactant = 1.0,
        enable_heat_transfer = True,
        enable_evaporation = True,
        enable_reactions = True,
        enable_dissipated_energy = True,
        excitation_parameters = [-2.0e5, 30000.0, 1.0],
        excitation_type = 'sin_impulse',
    )


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
    json_path: str,
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

    return dict(
        mechanism = 'chemkin_ar_he',
        R_E = dict(type="LinearRange", start=0.000005, end=0.000125, num_steps=30),
        species = ['O2'],
        fractions = [1.0],
        P_amb = dict(type="Const", value=101325.0),
        T_inf = dict(type="Const", value=293.15),
        alfa_M = dict(type="Const", value=0.35),
        P_v = dict(type="Const", value=2338.1),
        mu_L = dict(type="Const", value=0.001),
        rho_L = dict(type="Const", value=998.2),
        c_L = dict(type="Const", value=1483.0),
        surfactant = dict(type="Const", value=1.0),
        enable_heat_transfer = True,
        enable_evaporation = True,
        enable_reactions = True,
        enable_dissipated_energy = True,
        target_specie = 'H2',
        excitation_params = [
            dict(type="LinearRange", start=-150000.0, end=-250000.0, num_steps=20),
            dict(type="Const", value=30000.0),
            dict(type="Const", value=1.0)
        ],
        excitation_type = 'sin_impulse',
    )


def run_parameter_study(
    parameter_study: dict,
    json_path: str,
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
# User documentation **Bubble dynamics simulation C++**

This project is a high-performance sonochemistry simulation using C++20. The bubble dynamics, thermodynamics, and extremely high pressure and high temperature chemical kinetics of acoustically excited individual micron sized bubbles are simulated. The underlying ODE system is solved by the C-based [SUNDIALS CVODE](https://sundials.readthedocs.io/en/latest/cvode/Introduction_link.html) solver. This is a C++ powered version of the now obsolete [Bubble_dynamics_simulation](https://github.com/kozakaron/Bubble_dynamics_simulation).

To use this software, you have two options:
 1) **Build from source**: See [./CONTRIBUTING.md](./CONTRIBUTING.md) for details. This method is harder, but allows you more control.
 2) **Use the interface**: Download and unzip the latest release from GitHub Releases. You will find executables for both Windows (main.exe) and Linux (main). Use either the Python or MATLAB interface to run and plot individual simulations or to run and process brute-force parameter studies. Prerequisites for users: 
    * **Python 3**: For the interface. Remove affected functions from [./interface/python_interface.py](./interface/python_interface.py) to use without *pandas* and *matplotlib*.
        * **numpy**: Heavily used in the interface.
        * *pandas*: Used in the interface, but for reading back parameter studies only.
        * *matplotlib*: Used in the interface, but only in the `plot()` function.
    * *MATLAB*: Alternative for the interface.

 > [Video](https://youtu.be/XAMtw1O4HdY) for using the interface or building from source.

## Table of contents

- [Use the interface](#use-the-interface)
    - [Command line interface](#command-line-interface)
    - [Python interface](#python-interface)
    - [MATLAB interface](#matlab-interface)
- [Run individual simulations](#run-individual-simulations)
    - [example_cpar()](#example_cpar)
    - [run_simulation()](#run_simulation)
    - [plot()](#plot)
- [Run parameter studies](#run-parameter-studies)
    - [example_parameter_study()](#example_parameter_study)
    - [run_parameter_study()](#run_parameter_study)
    - [read_parameter_study()](#read_parameter_study)
    - [line_to_dict()](#line_to_dict)

## Use the interface

To use the interface, you simply have to download the latest release from GitHub. After unzipping it, you should find three subdirectories:
 * **bin**: Contains two executables: *main.exe* for Windows and *main* for Linux.
 * **mechanism**: Chemical mechanisms to simulate reactions for various bubble contents.
 * **interface**: Contains the MATLAB and Python interfaces.
 * You should place your MATLAB and Python files next to these subdirectories.

 > Note that all scripts use relative paths, so make sure that your current working directory is where this [./README.md](./README.md) is located. Your terminal should write: `some/path/Bubble_dynamics_simulation_CPP>`.

### Command line interface

 > You should use the Python or MATLAB based interfaces, skip this chapter.

You can use either the MATLAB or Python interfaces, however the latter one is recommended. Both of them are just wrappers for the command line interface. The interface utilizes the human readable JSON format. All settings are saved into a JSON file, such as [./interface/example_cpar.json](./interface/example_cpar.json) or [./interface/example_parameter_study.json](./interface/example_parameter_study.json). Then, the executable is called with command line arguments providing these JSON files.

The results of a single simulation are written back into the same JSON, adding to it. The numerical solution is dumped as unreadable binary data at the end of the JSON. Parsing time would be very long otherwise. An example command line command:
```
./bin/main.exe --run ./copy_of_example_cpar.json --tmax 1.0 --timeout 60.0 --save
```

Brute-force parameter studies are saved into their respective directories. The program will automatically number the directory to avoid name collisions. This directory will contain the simulation data in CSV files, the error log and the standard output, as well as the settings in both a JSON and a text file. An example command line command is:
```
./bin/main.exe --parameter_study ./interface/example_parameter_study.json --tmax 1.0 --timeout 60.0 --directory ./_parameter_studies/test
```

Use `./bin/main.exe --help` to refer to all options:
```
  -h, --help                 Show help
  -v, --version              Show version
  -r, --run arg              Run a simulation with the given JSON file
      --tmax arg             Simulation end time in seconds (default: 1.0)
      --timeout arg          Timeout in seconds (default: 60.0)
      --save                 Set this flag to save all timesteps, skip
                             it to save only the first and last steps
      --log arg              Set log file
      --parameter_study arg  Run a parameter study with the given JSON file
      --directory arg        Set save directory for parameter_study
                             (default: ./_parameter_studies/test)
```

### Python interface

You may import the interface and print the documentation with
```Python
from interface import python_interface as api
help(api)
# use as api.function_name()
```

Minimal Python example:
```Python
from interface import python_interface as api
cpar = api.example_cpar()
data = api.run_simulation(cpar)
api.plot(data)
```

### MATLAB interface

Use the following line to access the interface and print documentation:

```Matlab
addpath('interface');
help matlab_interface
% use as matlab_interface.function_name()
```

Minimal MATLAB example:
```Matlab
addpath('interface');
cpar = matlab_interface.example_cpar();
data = matlab_interface.run_simulation(cpar);
matlab_interface.plot(data);
```

## Run individual simulations

#### example_cpar()

Example codes use the Python interface. However, names are the same and functionalities are similar in MATLAB as well. To run and plot an individual parameter study, you should first assemble a dictionary containing the control parameters. To see all possible settings, **please refer to [./SIMULATION_SETTINGS.md](./SIMULATION_SETTINGS.md#inputs)**.

As a template, you may use `example_cpar()` and modify the returned dictionary:

```Python
cpar = api.example_cpar()
cpar['mechanism'] = 'chemkin_kaust2023_ammonia'
cpar['species'] = ['H2', 'N2']
cpar['fractions'] = [0.25, 0.75]
cpar['enable_evaporation'] = True
```

Similarly, in MATLAB:
```Matlab
cpar = matlab_interface.example_cpar();
cpar.mechanism = 'chemkin_kaust2023_ammonia';
cpar.species = {'H2', 'N2'};                   % species names in initial bubble
cpar.fractions = [0.25, 0.75];                 % molar fractions of species in initial bubble
cpar.enable_evaporation = true;                % toggles evaporation (ambient liquid is assumed to be water)
```

Or you may create a brand new dictionary from scratch:

```Python
cpar = {
    'ID': 0,
    'mechanism': 'chemkin_otomo2018_ammonia',
    'R_E': 1e-05,
    'ratio': 1.0,
    'species': ['H2', 'N2'],
    'fractions': [0.25, 0.75],
    'P_amb': 101325.0,
    'T_inf': 293.15,
    'alpha_M': 0.35,
    'P_v': 2338.1,
    'mu_L': 0.001,
    'rho_L': 998.2,
    'c_L': 1483.0,
    'surfactant': 1.0,
    'enable_heat_transfer': True,
    'enable_evaporation': True,
    'enable_reactions': True,
    'enable_dissipated_energy': True,
    'enable_van_der_waals': True,
    'bubble_dynamics': 'gilmore_nasg',
    'liquid_eos_params': [1.19, 6.218e8, 6.72e-4, 1.0e5, 997.0, 3610.0],
    'enable_rate_thresholding': True,
    'target_specie': 'H2',
    'excitation_type': 'sinusoid',
    'excitation_params': [-200000.0, 30000.0],
    'excitation_cycles': 1,
    'ramp_up_cycles': 0
 }
```

 > Incorrect types, invalid dictionary keys, unexpected values will trigger warnings or errors. Almost all mistakes are caught during preprocessing.

 > If you don't remember available chemical mechanisms, species, or excitation types/arguments, you can deliberately mistype  (e.g.: `cpar['mechanism'] = 'x'`) or set an invalid number of arguments (e.g.: `cpar['excitation_params'] = []`) and the resulting error message shall list the available options or the names/units of the arguments.


#### run_simulation()

In order to run the simulation, use `run_simulation()`:
```Python
data = api.run_simulation(cpar)
```
Arguments:
 * `cpar` (dict): Control parameters for the simulation. This is the only mandatory argument.
 * `json_path` (str): Path to the JSON file, which will be generated. (default: *ignore.json*)
 * `executable_path` (str): Path to the executable. (e.g.: *./bin/main.exe*)
 * `t_max` (float): Maximum simulation time. (default: 1.0 sec)
 * `timeout` (float): Timeout for the simulation. (default: 60.0 sec)
 * `save_steps` (bool): Whether to save the simulation steps. If false, only the first and last steps will be saved, which is a tiny bit faster. (default: true)
        
Returns (dict): A dictionary containing the simulation results. See key `'sol'` for the numerical solution and `'cpar'` for the original control parameters. The time series can be accessed via `data['sol']['t']` and `data['sol']['x']` as arrays. The simulation is successful if `data['sol']['success']` is true. You may also access post processing data under `data['postproc']`. Extra information is saved under keys `'excitation'`, `'mechanism'` and `'version'`.



#### plot()

You can create a plot about the temporal evolution of the bubble radius, internal temperature, and amount of compounds within the bubble. Also, information is printed in a more human readable manner. This function has somewhat less functionality in Matlab. Example usage of the `plot()` function:
```Python
api.plot(data)
```

Arguments:
 * `data` (dict): data dictionary (returned by `run_simulation()`). This is the only mandatory argument.
 * `n` (int): how long should the plotted time interval be compared to the collapse time. If the simulation is unsuccessful or n is negative, the entire simulation is plotted. (default: 5 [-])
 * base_name: save plots as images (default: '' alias do not save).
    * use base_name='plot' --> plot_1.png, plot_2.png 
    * use base_name='images/plot' to save into images folder. Using a folder for images is recommended. This folder has to be created manually
 * `format`: format of the saved images (available: png, pdf, ps, eps, svg)
 * `presentation_mode`: if True, the plot will be in presentation mode, with thicker lines and larger fonts. (default: False)
 * `show_legend`: if True, the legend will be visible with every single species (default: False)
 * `show_cpar`: if True, the control parameters will be printed on the plot (default: False)
 * `plot_pressure`: if True, the excitation and internal pressure will be plotted as well (default: False)

## Run parameter studies

A brute-force parameter study automatically runs several simulations in parallel, and saves only the post processing data in CSV format.

#### example_parameter_study()

 Similarly to individual simulations, a parameter study is initiated by a dictionary. You may get a template with `example_parameter_study()`:
 ```Python
parameter_study = api.example_parameter_study()
 ```

Just like the case of the control parameters, you need to provide the mechanism and excitation type as strings, enable_ variables as bool and species and fractions as lists. However, numerical variables are provided as a range, see [./include/parameter_study.h](./include/parameter_study.h). There are 3 available options:
 * **Constant**: This variable is set as a fixed value, and it isn't changed in the parameter study. E.g.: `{"type": "Const", value=1.0}`
 * **Linear range**: This parameter takes part in the parameter study. The variable is changed from start to stop with num_step equal increments. E.g.: `{"type": "LinearRange", start=0.0, end=1.0, num_steps=10}`
 * **Logarithmic range**: This parameter takes part in the parameter study. The variable is changed from start to stop with num_step increments. Subdivision is uneven, changing the same way as Python's numpy.logspace. GeomRange is more versatile. E.g.: `{"type": "LogRange", start=1.0, end=100.0, num_steps=10}`
 * **Geometric series range**: This parameter takes part in the parameter study. The variable is changed from start to stop with num_step increments. Subdivision is uneven, the difference of consecutive elements forms a geometric series with a quotient q. If 0 < q < 1, steps are getting gradually smaller. If 1 < q, steps are getting larger. Increase q to get a larger increase in steps. E.g.: `{"type": "GeomRange", start=1.0, end=100.0, num_steps=10, q=2.0}`

 A sample of such a control dict may look like this:
 ```Python
 parameter_studies = {
    'mechanism': 'chemkin_elte2016_hydrogen',
    'R_E': {
        'type': 'LinearRange',
        'start': 5e-06,
        'end': 0.000125,
        'num_steps': 30
    },
    'species': ['O2'],
    'fractions': [1.0],
    'P_amb': {'type': 'Const', 'value': 101325.0},
    # ...
    'excitation_params': [
        {
            'type': 'LinearRange',
            'start': -150000.0,
            'end': -250000.0,
            'num_steps': 20
        },
        {'type': 'Const', 'value': 30000.0},
        {'type': 'Const', 'value': 1.0}],
    'excitation_type': 'sinusoid'
}
 ```

#### run_parameter_study()

Once you have your settings, you can call the parameter study with `run_parameter_study()`:
```Python
api.run_parameter_study(parameter_study)
```

Arguments:
 * `parameter_study` (dict): Parameters for the simulation.
 * `json_path` (str): Path to the JSON file. (default: *ignore.json*)
 * `executable_path` (str): Path to the executable. (e.g.: *./bin/main.exe*)
 * `t_max` (float): Maximum simulation time. (default: 1.0 sec)
 * `timeout` (float): Timeout for the simulation. (default: 60.0 sec)
 * `save_directory` (str): Where to save the parameter study. (default: *./_parameter_studies/test*)

 A parameter study will automatically create a directory at `save_directory`. If the path already exists, it will automatically number the parameter studies, e.g.: *test1*, *test2*, ... . This directory shall contain one CSV file for each thread containing the simulation results in a tabular manner. Also, the following files:
```
_parameter_studies
 ├─test
 │ └─... // older study
 └─test2
    ├─sundials_logs
    | └─... // possible SUNDIALS warning, which may clutter console output if printed there.
    ├─bruteforce_parameter_study_settings.json  // json to reproduce the parameter study with the interface + info
    ├─bruteforce_parameter_study_settings.txt   // same, but as C++ code
    ├─errors.log    // a log file listing all errors/warnings, which occurred during the parameter study
    ├─output_0.csv  // results of thread 0
    ├─output_1.cev  // results of thread 1
    ├─...
    ├─output.log    // a log file of the standard output
    └─summary.log   // summary of the parameter study, if it finished (runtime, number of failed simulations)
```

#### read_parameter_study()

You may read the results of a parameter study into a pandas dataframe with `read_parameter_study()`. This dataframe is sorted according to `energy_demand`:
```Python
all_data = read_parameter_study('path/to/parameter/study')
good_data = all_data.loc[all_data['success'] == True]   # locate successful simulations
```

#### line_to_dict()

To investigate a single simulation from a parameter study, use `line_to_dict()`:
```Python
cpar = api.line_to_dict(all_data.iloc[0])   # best energy demand
data = api.run_simulation(cpar)
api.plot(data)
```

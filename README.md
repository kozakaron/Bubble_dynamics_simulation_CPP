# Bubble dynamics simulation C++
A C++ powered version of [Bubble_dynamics_simulation](https://github.com/kozakaron/Bubble_dynamics_simulation).

To use this software, you have two options:
 1) **Build from source**: Clone the git repository, run the [./dev/build_submodules.py](./dev/build_submodules.py) script to check dependencies and install external repositories. Then you may run [./dev/build.py](./dev/build.py) to compile the code. This method is harder, but allows you more control.
 2) **Use the inteface**: Download and unzip the latest release from GitHub Releases. You will find executables for both Windows (main.exe) and Linux (main). Use either the Python or Matlab interface to run and plot indavidual simulations or to tun and process bruteforce parameter studies.

Note, that all scripts use relative paths, so make sure that your current working directory is where this [./README.md](./README.md) is located. Your terminal should write: `some/path/Bubble_dynamics_simulation_CPP>`.

## Build from source

### Cloning repository

You choosed the harder path. To build from source, you must download this repository via git, and not as zip or trough any other means. Use:
```
git clone https://github.com/kozakaron/Bubble_dynamics_simulation_CPP.git
```

### Prerequisits

The first step is to check prerequisits. There is a costum script, [./dev/build_submodules.py](./dev/build_submodules.py), you should run:
```
python ./dev/build_submodules.py
```
This should scan requirements, and print something like `git: ok; cmake: ok; gcc: ok; ...`. You should install missing softwares, either by the provided link or command. You only need one compiler (gnu g++ or llvm clang++) and the debuggers (gdb or lldb) are optional. For the build process, you need Visual Studio on Windows and make on Linux. All other requirements are mandatory.

### Installing submodules

After the requirements, the submodules should be cloned and installed. These are third party GitHub repositories. I used GitHub's submodules functionality to link them into my repository. The same script clones these external repositories, and builds them. When the scripts asks `Do you want to build SUNDIALS from source? [y/n]`, enter y. This may take a few minutes, and problems can arise. You may modify build setting at the top of [./dev/build_submodules.py](./dev/build_submodules.py) by changing `cmake_command_list` in the `submodules` dictionary. You may also have to modify `vs_path` on Windows by providing your Visual Studio installation path. Make sure, that neither submodule failed to install, otherwise the linking process will fail.

### Compilation

 I hate make and cmake, thus a python script is used to compile the project.
 You may change the include directories, compiler, debugger, and related flags in the [./dev/build.py](./dev/build.py) script. However, beside turning off the sanitizers for high performance benchmarking, you will probably just use the script's flags:
~~~
  -h, --help               show this help message and exit
  -t, --test               build with test flag (TEST is defined)
  -b, --benchmark          build with benchmark flag (BENCHMARK is defined)
  -r, --run                run the binary after building
  -w, --warning            enable all warnings
  -d, --debug              compile in debug mode (also runs gdb if called with -r)
  -s, --shared             build as shared library (.dll/.so)
  -o2, -O2, --optimize2    enable o2 optimization
  -o3, -O3, --optimize3    enable o3 optimization
~~~
 The script collects each *.cpp* file automatically, and compiles them into object files, each one on a seperate thread. Then the object files are linked togeather along with the external repositories' lib files. Syntax colored error messages are printed grouped by *.cpp* file. Use flags `--test`, `--benchmark` to run unit tests or benchmarks respectively. Use `--run` to run the executable automatically after compilation. 

For example you may compile+run your tests and benchmarks as: `python ./dev/build.py --run --warning --optimize3 --benchmark --test` or `python ./dev/build.py -r -w -o3 -b -t`. Check the [./bin](./bin) directry for the executables and the build logs. The latter one should contain all errors, the commands used to compile each object file and for linking, also the build times.

Notes:
 * The main entry point is in [./src/main.cpp](./src/main.cpp).
 * You may have problem with the **sanitizers**, especially on Windows. Feel free to turn them off. Also turn them off for high performance computing. However, they are very powerful tools to find bugs. I never managed to make the thread sanitizer work.
 * Some systems strugle with the long double. If the long double exists unit test fails or you encounter errors, you may try to change `-mlong-double-80` to 64 or 128. (I actually found a compiler error, causing segfault while calculating with long doubles on windows: [link](https://github.com/llvm/llvm-project/issues/118138))
 * Thanks to multithreading, compilation performance is very good. However, in the future, it can be further improved. First, by using precompiled headers. And second, by creating a hash of each preprocessed source file, and only recompileing the ones that changed.

 ## Use the interface

To use the interface, you simply have to download the latest release from GitHub. After unzipping it, you should find two directories:
 * **bin**: contains two executables: *main.exe* for Windows and *main* for Linux.
 * **interface**: contains the Matlab and Python interfaces.

### Command line interface

You can use either the Matlab or Python interfaces, however the latter is more recommended. Both of them are just wrappers for the command line interface. The interface utilize the human readable JSON format. All settings are saved into a JSON file, such as [./interface/example_cpar.json](./interface/example_cpar.json) or [./interface/example_parameter_study.json](./interface/example_parameter_study.json). Then, the executable is called with command line arguments providing these JSON file.

The results of a single simulation are written back into the same JSON, adding to it. The numerical solution is dumped as unreadable binary data at the end of the JSON. Parsing time would be very long otherwise. An example command line command:
```
./bin/main.exe --run ./copy_of_example_cpar.json --tmax 1.0 --timeout 60.0 --save
```

Bruteforce parameter studies are saved into their respective directories. The program will automatically number the directory to avoid name collisions. This directory will contain the simulation data in csv files, the error log and the standard output, as well as the settings in both a json and a text file. An example command line command is:
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

 ### Matlab interface

 Use the following line to access the interface and print documentation:

 ```Matlab
 addpath('interface');
 help matlab_interface
 % use as matlab_interface.function_name()
 ```

 ### Run indavidual simulations

 #### example_cpar()

Example codes use the Python interface. However, names are the same and functionalities are similar in Matlab as well. To run and plot an indavidual parameter study, you should first assemble a dictionary containing the control parameters. As a template, you may use `example_cpar()` and modify the returned dictionary:
```Python
cpar = api.example_cpar()
cpar['mechanism'] = 'chemkin_otomo2018'
cpar['species'] = ['H2', 'N2']
cpar['fractions'] = [0.25, 0.75]
cpar['enable_evaporation'] = True
```

Or you may create a brand new dictionary:
```Python
cpar = {
    'ID': 0,
    'mechanism': 'chemkin_otomo2018',
    'R_E': 1e-05,
    'species': ['H2', 'N2'],
    'fractions': [0.25, 0.75],
    'P_amb': 101325.0,
    'T_inf': 293.15,
    'alfa_M': 0.35,
    'P_v': 2338.1,
    'mu_L': 0.001,
    'rho_L': 998.2,
    'c_L': 1483.0,
    'surfactant': 1.0,
    'enable_heat_transfer': True,
    'enable_evaporation': True,
    'enable_reactions': True,
    'enable_dissipated_energy': True,
    'target_specie': 'H2',
    'excitation_params': [-200000.0, 30000.0, 1.0],
    'excitation_type': 'sin_impulse'
 }
```

If you are unsure about the available mechanisms, excitation types or excitation parameters, you may refer to [./include/parameters.h](./include/parameters.h). As an alternative, you can mistype the type (e.g.: `cpar['mechanism'] = 'x'`) or set an invalid number of arguments (e.g.: `cpar['excitation_params'] = []`) and the resulting error message shall list the available options or the names/units of the arguments.

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
 * `save_steps` (bool): Whether to save the simulation steps. If false, only the first and last steps will be saved, which is a tony bit faster. (default: true)
        
Returns (dict): A dictionary containing the simulation results. See key `'sol'` for the numerical solution and 'cpar' for the original control parameters. the numerical solution can be accessed via `data['sol']['t']` and `data['sol']['x']` as arrays. The simulation is sucessful if `data['sol']['success']` is true. You may also acces post processing data, like `data['energy_demand']`. Extra information is saved under keys `'excitation'`, `'mechanism'` and `'version'`.

#### plot()

You can create a plot about the temporal evolution of the bubble radius, internal temperature, and amount of compounds within the bubble. Also, information is printed in a more human readable manner. Use the `plot` function. This function has somewhat less functionality in Matlab. Example usage:
```Python
api.plot(data)
```

Arguments:
 * `data` (dict): data dictionary (returned by `run_simulation()`). This is the only mandatory argument.
 * `n` (int): how long should the plotted time interval be compared to the collapse time. If the simulation is unsuccessful or n is negative, the entire simulation is plotted. (default: 5 [-])
 * base_name: save plots as image (default: '' alias do not save).
    * use base_name='plot' --> plot_1.png, plot_2.png 
    * use base_name='images/plot' to save into images folder. Using a folder for images is recommend. This folder have to be created manually
 * `format`: format of the saved images (available: png, pdf, ps, eps, svg)
 * `presentation_mode`: if True, the plot will be in presentation mode, whith thicker lines and larger fonts. (default: False)
 * `show_legend`: if True, the legend will be visible with every single species (default: False)
 * `show_cpar`: if True, the control parameters will be printed on the plot (default: False)

 ### Run parameter studies

A bruteforce parameter study automatically runs several simulations in parallel, and saves only the post processing data in CSV format. When you iniciate a  

#### example_parameter_study()

 Similarly to indavidual simulations, a parameter study is initiated by a dictionary. you may get a template with `example_parameter_study()`:
 ```Python
parameter_study = api.example_parameter_study()
 ```

Just like as for control parameters, you need to provide the mechanism and excitation type as strings, enable_ variables as bool and species and fractions as lists. However, numerical variables are provided as a range, see [./include/parameter_study.h](./include/parameter_study.h). There are 3 available options:
 * **Constant**: This variable is set as a fixed value, isn't changed in the parameter study. E.g.: `{"type": "Const", value=1.0}`
 * **Linear range**: This parameter takes part in the parameter study. The variable is changed from start to stop with num_step equal increment. E.g.: `{"type": "LinearRange", start=0.0, end=1.0, num_steps=10}`
 * **Power range**: This parameter takes part in the parameter study. The variable is changed from start to stop with num_step increments. Subdivision is uneven, controlled by power. E.g.: `{"type": "PowRange", start=1.0, end=100.0, num_steps=10, base=2.0}`

 A sample of such a control dict may look like this:
 ```Python
 parameter_studies = {
    'mechanism': 'chemkin_ar_he',
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
    'excitation_type': 'sin_impulse'
}
 ```

#### run_parameter_study()

Once you have your settings, you can call the parameter study with `run_parameter_study()`:
```Python
api.run_parameter_study(
    parameter_study,
    json_path='ignore.json'
)
```

Arguments:
 * `parameter_study` (dict): Parameters for the simulation.
 * `json_path` (str): Path to the JSON file. (default: *ignore.json*)
 * `executable_path` (str): Path to the executable. (e.g.: *./bin/main.exe*)
 * `t_max` (float): Maximum simulation time. (default: 1.0 sec)
 * `timeout` (float): Timeout for the simulation. (default: 60.0 sec)
 * `save_directory` (str): Where to save the parameter study. (default: *./_parameter_studies/test*)

 A parameter study will automatically create a directory at `save_directory()`. If the path already exists, it will automatically number the parameter studies, e.g.: *test1*, *test2*, ... . This directory shall contain one CSV file for each thread containing the simulation results in a tabular manner. Also, the following files:
```json
_parameter_studies
 ├─test
 │ └─... // older studies
 └─test2
    ├─bruteforce_parameter_study_settings.json  // json to reproduce the parameter study with the interface + info
    ├─bruteforce_parameter_study_settings.txt   // same, but as C++ code
    ├─errors.log    // a log file listing all errors/warnings, which occured during the parameter study
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

## How it works?

This code is object oriented. Most classes include some of the following methods. `ControlParameters` from [./include/control_parameters.h](./include/control_parameters.h) is an excellent example.
 * **constructor with a builder**: many classes have a builder struct, like `ControlParameters::Builder`. This struct allows the usage of designated initializer lists. You no longer have to know the order of constructor arguments, and default values are handled as well.
 * **constructor from JSON**: some classes can be constructed from a JSON to use with the interface. They accept both the file paths for the JSON and `nlohmann::ordered_json`. Some class also has a `to_json()` method to help save as file.
 * **to_string()**: use this to convert the object to a human readable string. Some arguments might be available, such as `colored` to print with some syntax coloring. Or `with_code`, which helps to print objects as valid C++ code. This is done in the form of builder struct, which you can copy into your own code.
 * **operator<< overload**: use to print objects directly trough `std::cout`.
 * **to_csv()**: convert the object to a csv (sting). Member variables are printed, seperated by a comma. Order is defined by `csv_header`.

### Error handling

See `Error` and `ErrorHandler` classes in [./include/common.h](./include/common.h). Each error has:
 * A severity (info, warning, error) and a type (general, preprocessing, ...). The type is meant to be used for statistical handling or simulation failiures.
 * A string error meassage, usually with detailed context.
 * Further metadata: exact timedate, location (file, function, line), and ID of the related control parameters combination, if the error is related to a specific simulation.

Errors are logged with the `LOG_ERROR(...)` macro, e.g.: `size_t error_ID = LOG_ERROR(Error::severity::error, Error::type::preprocess, message_string);`. Some classes have a member called `error_ID` used to store and propagate a state of error. `ErrorHandler::no_error` is used if the object is functioning as expected. Errors are printed as:

```
2025.02.15 21:31:08: (general) here is the error message with details ./src/some_file.cpp:42: Class::function()
```

`ErrorHandler` is thread safe. Use `ErrorHandler::set_log_file()` to log into a file, and set `ErrorHandler::print_when_log = false;` to disable printing on the console.

### Chemical mechanisms

Chemical mechanisms are generated from chemkin *.inp* files. (Not included, see [list of .inp files](./dev/inp_file_list.py)) They are turned into usable format via the (quite messy) [./dev/inp_data_extractor.py](./dev/inp_data_extractor.py). The generated structs in [./mechanism](./mechanism) have static constexpr members only. In order to swithc between mechanisms in runtime, and make migrating to CUDA easier, these structs are turned into const static members of the `Parameters` class in [./include/parameters.h](./include/parameters.h). This class stores all arrays in a flattened form as raw pointers. Usage: `const Parameters *par = Parameters::get_parameters(Parameters::mechanism::chemkin_ar_he);` and get any members like `par->nu`.

### Control parameters, and ODE function

The `ControlParameters` class, declared in [./include/control_parameters.h](./include/control_parameters.h), holds all parameters, which may influence the simulation. The recommended way to initialize is with builder struct and a designated initializer list: 

```cpp
ControlParameters cpar = ControlParameters{ControlParameters::Builder{
    .ID                          = 0,
    .mechanism                   = Parameters::mechanism::chemkin_ar_he,
    .R_E                         = 1.00000000000000008e-05,    // bubble equilibrium radius [m]
    .species                     = {"O2"},
    .fractions                   = {1.00000000000000000e+00},
    .P_amb                       = 1.01325000000000000e+05,    // ambient pressure [Pa]
    .T_inf                       = 2.93149999999999977e+02,    // ambient temperature [K]
    .alfa_M                      = 3.49999999999999978e-01,    // water accommodation coefficient [-]
    .P_v                         = 2.33809999999999991e+03,    // vapour pressure [Pa]
    .mu_L                        = 1.00000000000000002e-03,    // dynamic viscosity [Pa*s]
    .rho_L                       = 9.98200000000000045e+02,    // liquid density [kg/m^3]
    .c_L                         = 1.48300000000000000e+03,    // sound speed [m/s]
    .surfactant                  = 1.00000000000000000e+00,    // surface tension modifier [-]
    .enable_heat_transfer        = true,
    .enable_evaporation          = true,
    .enable_reactions            = true,
    .enable_dissipated_energy    = true,
    .target_specie               = "H2",
    .excitation_params           = {-2.00000000000000000e+05, 3.00000000000000000e+04, 1.00000000000000000e+00},
    .excitation_type             = Parameters::excitation::sin_impulse
}};
```
These above are also the default values, any builder argument can be missed. Use the `to_sting()` method or the `ostream operator<<` overload to print the control parameters to the console. It will be printed in the same format as above, which is also valid code. The class also have a csv header and a `to_csv()` method, like many other classes.

The right hand side of the ODE is calculated by the `OdeFun` class declared in [./include/ode_fun.h](./include/ode_fun.h). (The right hand side if $f()$ if the ODE is $\frac{dx}{dt} = f(x)$ and the initial condition is $x(0) = x_0$) `OdeFun` is a function object, with the call operator (`operator()`) defined. It has to be initialized, so it can allocate and reuse arrays for temporary results, and not call `new` during high performance calculations:

```cpp
OdeFun ode;
ode.init(cpar);
```

`OdeFun` is also responsible to calculate the initial conditions, however, an array of appropiate size has to be provided as a pointer. This step is handled by the solver:

```cpp
std::vector<double> x_0(ode.par->num_species+4);
ode.initial_conditions(x_0.data());
```

Now, the right hand side can be calculated with `operator()`. The pointers for $x$ and $\frac{dx}{dt}$ has to be provided. This function returns the type `is_success`, which is effectievly a boolean. However, in case of an error, `ode.cpar.error_ID` will also change from `Error::no_error`. You may obtain the error as:
```cpp
Error error = ErrorHandler.get_error(ode.cpar.error_ID);
std::cout << error <<std::endl;
```

### Runing a simulation

Note: This part is expected to change.

```cpp
// init cpar and ode
ControlParameters cpar = ControlParameters{{ /* ... */ }};
OdeFun ode;
ode.init(cpar);

// solve the ODE      
OdeSolverCVODE solver(ode.par->num_species+4);
OdeSolution solution = solver.solve(
    1.0,     // t_max [s]
    &ode,    // ode_ptr
    60.0,    // timeout [s]
    true     // weither to save solution (or just first and last step)
);

// do postprocessing and print results (with operator<< overload)
SimulationData data(cpar, sol);
std::cout << data << std::endl;
```

### Running a bruteforce parameter study

A bruteforce parameter study can be defined by the `ParameterCombinator` class declared in [./include/parameter_study.h](./include/parameter_study.h). It is also initialized with a builder struct and a designated initializer list, similar to `ControlParameters`. However, arguments have to be a childre of the `Range` class: `Const`, `LinearRange`, `PowRange`. It also has defaults for missing arguments, and can be printed to console in usable code format:

```cpp
ParameterCombinator parameter_combinator = ParameterCombinator{ParameterCombinator::Builder{
    .mechanism                   = Parameters::mechanism::chemkin_ar_he,
    .R_E                         = LinearRange(0.000005, 0.000125, 5),                 // {5e-06, 3.5e-05, 6.5e-05, 9.5e-05, 0.000125}
    .P_amb                       = Const(101325.000000),                               // {101325}
     /* ... */
    .enable_dissipated_energy    = true,
    .target_specie               = "H2",
    .excitation_params           = {
        LinearRange(-100000.000000, -300000.000000, 5),     // {-100000, -150000, -200000, -250000, -300000}
        Const(20000.000000),                                // {20000}
        Const(1.000000)                                     // {1}
    },
    .excitation_type             = Parameters::excitation::sin_impulse
}};

std::cout << parameter_combinator << std::endl;
```

The total number of all possible control parameter combinations are given by `parameter_study.get_total_combination_count()`. You may iterate trough these combinations in a thread safe way as:

```cpp
while (true)
{
   auto [success, cpar] = parameter_study.get_next_combination();
   if (!success) break;
   /* solve ODE corresponding to cpar */
}
```

The `ParameterCombinator` class does not create all possible `ControlParameters` ahead of time in a huge list. Thus, it can theoritically handle an arbitrary number of combinations with a small memory footprint.

Use `ParameterStudy` class to create save folder, run simulations multithreaded, log output and errors, save results to csv automatically. See [./test/benchmark_parameter_study.cpp](./test/benchmark_parameter_study.cpp) for example code.
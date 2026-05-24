# Developer guide for **Bubble dynamics simulation C++**

> Provide this file for AI coding agent, when used.

This project is a high-performance sonochemistry simulation using C++20. The bubble dynamics, thermodynamics, and extremely high pressure and high temperature chemical kinetics of acoustically excited individual micron sized bubbles are simulated. The underlying ODE system is solved by the C-based [SUNDIALS CVODE](https://sundials.readthedocs.io/en/latest/cvode/Introduction_link.html) solver. This is a C++ powered version of the now obsolete [Bubble_dynamics_simulation](https://github.com/kozakaron/Bubble_dynamics_simulation).

## Table of contents

- [Coding guidelines](#coding-guidelines)
- [Directory layout](#directory-layout)
- [Build from source](#build-from-source)
    - [Prerequisites for developers](#prerequisites-for-developers)
    - [Build prerequisites (first time build)](#build-prerequisites-first-time-build)
    - [Build the project](#build-the-project)
    - [Make a release, commit](#make-a-release-commit)
- [How it works?](#how-it-works)
    - [High level overview](#high-level-overview)
        - [Running a single simulation](#running-a-single-simulation)
        - [Running a parameter study](#running-a-parameter-study)
    - [Error handling](#error-handling)
    - [Control parameters and ODE function](#control-parameters-and-ode-function)
    - [Running a simulation](#running-a-simulation)
    - [Running a brute-force parameter study](#running-a-brute-force-parameter-study)
    - [Testing and benchmarking with built in unit testing framework](#testing-and-benchmarking-with-built-in-unit-testing-framework)

## Coding guidelines

 * Naming conventions: `ClassesAndStructs`, `functions_and_variables`. Use explanatory naming, group variables with namespaces and structs. E.g. when naming velocity components, instead of `u`, `v`, `w`, use `vel.x`, `vel.y`, `vel.z`.
 * Types: use types from `stdint.h`, e.g. `uint64_t` instead of `unsigned long long`. Use `bool` from `stdbool.h`.
 * Constants: Try to use `constexpr` over `#define`. Use `numbers.h`, e.g. `std::numbers::pi` instead of `MATH_PI`. When adding physical constants (like universal gas constant), store it in class `Parameters`.
 * Use object-oriented programming.
 * In high performance code ([./src/ode_fun.cpp](./src/ode_fun.cpp)), use raw pointers.
    * Resource allocation and freeing should be tied to variable lifespan using class constructors and destructors. (RAII)
    * Multidimensional arrays are stored in a flattened form. Instead of `array[i][j]` use `array[i*size_x+j]`
 * Otherwise use std containers (like `std::vector`) and smart pointers (`std::unique_ptr`, `std::shared_ptr`).
 * Error handling: don't use `try{...} catch(...){...}`, use the `LOG_ERROR` macro and store a state of error in the `error_ID` member.
    * Instead of `throw std::runtime_error("error");`, use `size_t error_ID = LOG_ERROR("error");`. You may ignore the returned value.
    * Or you may use `some_object.error_ID` to store a state of error. The value during normal functioning is `ErrorHandler::no_error`.
    * Functions or methods may return an `is_success` (bool) type, which is true during normal functioning, and false if an error occurred. In this case, set `error_ID`.
 * You build with `python ./dev/build.py -r -w -o2`. When developing from Linux or WSL, you may use address sanitizers, enable them in [./dev/build.py](./dev/build.py) under `compiler_flags`.
 * Try not to use system-specific code. No `windows.h`, no `pthread.h`, instead use `std::thread` (from `<thread>`).
 * Test everything! Even though unit tests might not have been updated :)
 * Keep formulas readable. Never trust AI to do math correctly. An example: $p_{ex}(t) = P_{\infty} + P_{A1} \cdot sin(2 \pi f_1 t) +  P_{A2} \cdot sin(2 \pi f_2 t + \theta) $
    ```C++
    // unreadable: (but short)
    const double p_excitation = cpar.P_amb + cpar.excitation_params[0] * sin(2.0 * MATH_PI * cpar.excitation_params[2] * t) + cpar.excitation_params[1] * sin(2.0 * MATH_PI * cpar.excitation_params[3] * t + cpar.excitation_params[4]);

    // readable:
    const double& p_A1        = cpar.excitation_params[0];
    const double& p_A2        = cpar.excitation_params[1];
    const double& freq1       = cpar.excitation_params[2];
    const double& freq2       = cpar.excitation_params[3];
    const double& theta_phase = cpar.excitation_params[4];

    const double p_ex1        = p_A1 * sin(2.0 * std::numbers::pi * freq1 * t);
    const double p_ex2        = p_A2 * sin(2.0 * std::numbers::pi * freq2 * t + theta_phase);
    const double p_excitation = cpar.P_amb + p_ex1 + p_ex2;
    ```

## Directory layout

Here is a brief overview of the project's structure to help you navigate:

| Directory | Description |
|:---|:---|
| `bin/` | Created dynamically. Output location for compiled executables, object files, and build logs (`build.log`). |
| `dev/` | Developer scripts and assets. Contains Python build scripts (`build.py`, `build_submodules.py`), and documentation images. |
| `include/` | Contains all the C++ headers (`.h`). The definitions for classes like `ControlParameters`, `OdeSolution`, and `OdeSolver`. |
| `interface/` | The MATLAB and Python interfaces used to call the simulation from higher level languages. Also includes example JSON files. |
| `mechanism/` | Folder containing chemical mechanisms. Stores the original `.inp` (Chemkin) files, converted `.yaml` (Cantera), in-project `.json` representations. Has scripts to help add new mechanisms. |
| `src/` | Contains all the C++ source files (`.cpp`) that implement the physics, solver interface, logging, and parameter handling. |
| `submodules/` | External third-party libraries (installed via Git submodules). Contains SUNDIALS CVODE, nlohmann JSON, and cxxopts. |
| `test/` | Contains C++ files that define unit tests, speed benchmarks. |

Key source files:

* `main.cpp`: The main entry point of the executable. Handles command-line arguments, parameter initialization, and calls the simulation.
* `common.h`, `common.cpp`: Contains the `ErrorHandler` class, logging utilities, and global helper functions.
* `ode_fun.h`, `ode_fun.cpp`: This is where you will find all the formulas implemented.
* `ode_solver.h`, `ode_solver.cpp`: The wrapper around the SUNDIALS CVODE library. Manages the solver setup, integration steps, and tolerances.
* `parameters.h`, `parameters.cpp`: This is where all the constants are stored, after loaded from [./mechanism/json_files/](./mechanism/json_files/).

## Build from source

### Prerequisites for developers

**mandatory**, *optional*

 * **Git**: To clone the repository and fetch its submodules.
 * **Python 3**: For the interface, and to run the custom `build.py` and `build_submodules.py` scripts.
    * **pygments**: Used in [./dev/build_utility.py](./dev/build_utility.py) to syntax-highlight the compiler's C++ error logs in the terminal.
    * **numpy**: Heavily used in the interface.
    * *pandas*: Used in the interface, but for rading back parameter studies only.
    * *matplotlib*: Used in the interface, but only in the `plot()` function.
 * *MATLAB*: Alternative for the interface.
 * **C++20 Compiler**: Choose one to build the project, should be specified in variable `cpp_compiler` in [./dev/build.py](./dev/build.py)
    * gnu g++: Great option, can be installed on Windows with the MinGW-w64 project.
    * llvm clang++: Great option, recommended for Linux.
 * *C++ debugger*: Specify in variable `debugger` in [./dev/build.py](./dev/build.py) to automatically launch with command line debugger, when the `--run`, `--debug` flags are provided.
    * *gnu gdb*: For g++.
    * *llvm lldb*: For clang++.
 * **CMake**: To build CVODE.
 * **Build tools for CMake**:
    * *make*: For Linux builds.
    * *Visual Studio*: CMAKE generates VS project files under Windows. Specify VS install path in variable `vs_path` in [./dev/build_submodules.py](./dev/build_submodules.py).
 * **C++ Submodules**: Automatically downloaded from GitHub by [./dev/build_submodules](./dev/build_submodules.py)
    * SUNDIALS CVODE: High performance and modular ODE solver. Can be tricky to build.
    * nlohmann's json: Header only library to parse JSON files.
    * jarro2783's cxxopts: Header only library to parse command line arguments.
    
> It is recommended to use Windows Subsystem for Linux (WSL) or Chocolatey Package Manager on Windows to install prerequisites. 

### Build prerequisites (first time build)

 1. To build from source, you must download this repository via git, and not as zip or through any other means. Use:
    ```
    git clone https://github.com/kozakaron/Bubble_dynamics_simulation_CPP.git
    ```
 2. The first step is to check prerequisites. There is a custom script, [./dev/build_submodules.py](./dev/build_submodules.py), which you should run:
    ```
    python ./dev/build_submodules.py
    ```
    This should scan requirements, and print something like `git: ok; cmake: ok; gcc: ok; ...`. You should install missing software, either by the provided link or command. You only need one compiler (gnu g++ or llvm clang++) and the debuggers (gdb or lldb) are optional. For the build process, you need Visual Studio on Windows and make on Linux. All other requirements are mandatory.
 3. After the requirements, the submodules should be cloned and installed. These are third party GitHub repositories. I used GitHub's submodules functionality to link them into my repository. The same script clones these external repositories, and builds them. When the script asks `Do you want to build SUNDIALS from source? [y/n]`, enter y. This may take a few minutes, and problems can arise. You may modify build settings at the top of [./dev/build_submodules.py](./dev/build_submodules.py) by changing `cmake_command_list` in the `submodules` dictionary. You may also have to modify `vs_path` on Windows by providing your Visual Studio installation path. Make sure that neither submodule failed to install, otherwise the linking process will fail later.
 4. Build the project and run with: (see next chapter for details)
    ```
    python ./dev/build.py -r -w -o2
    ```
    

### Build the project

 I hate make and cmake, thus a python script is used to compile the project. Both Windows and Linux are supported. 
 You may change the include directories, compiler, debugger, and related flags in [./dev/build.py](./dev/build.py). However, beside turning off the sanitizers for high performance benchmarking, you will probably just use the script's flags:
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
 The script collects each *.cpp* file automatically, and compiles them into object files, each one on a separate thread. Then the object files are linked together along with the external repositories' lib files. Syntax colored error messages are printed grouped by *.cpp* file. Use flags `--test`, `--benchmark` to run unit tests or benchmarks respectively. Use `--run` to run the executable automatically after compilation. 

For example you may compile+run your benchmarks as: `python ./dev/build.py --run --warning --optimize3 --benchmark` or `python ./dev/build.py -r -w -o3 -b`. Check the [./bin](./bin) directory for the executables and the build logs (*./bin/build.log*). The latter one should contain all errors, the commands used to compile each object file and for linking, also the build times.

Notes:
 * All scripts use relative paths, so make sure that your current working directory is where this [./README.md](./README.md) is located. Your terminal should write: `some/path/Bubble_dynamics_simulation_CPP>`.
 * The main entry point is in [./src/main.cpp](./src/main.cpp).
 * You may have problems with the **sanitizers**, especially on Windows. Feel free to turn them off. Also turn them off for high performance computing. However, they are very powerful tools to find bugs. I never managed to make the thread sanitizer work.
 * Some systems struggle with the long double. If the long double exists unit test fails or you encounter errors, you may try to change `-mlong-double-80` to 64 or 128.
 * Unit tests won't build in the current version.


### Make a release, commit

When working on the source code, work on your own branch, the main branch should have a working version. Have someone review your pull requests. Please name your branches, and write a brief commit message. Please only commit files, which are absolutely needed for the core project, do not commit test Python scripts, simulation results, logs, excel files, large files. Use the [./.gitignore](./.gitignore).

To make a release on GitHub, make a new directory with the correct version number. For non-beta releases, increment the version number in [./include/common.h](./include/common.h). Copy directories *interface* (without *\_\_pycache__*), and *mechanism*, markdown files, and make an empty *bin* subdirectory. Build the submodules and the executables (with O2 optimization) on both Windows and Linux (WSL). Copy just the executables to the *bin* subdirectory. Zip the directory and upload on GitHub under releases. 

```
Bubble_dynamics_simulation_CPP_x.x  # replace x with version number, e.g.: Bubble_dynamics_simulation_CPP_1.7
├── bin                             # delete logs and object files
│   ├── main
│   └── main.exe
├── interface                       # Delete __pycache__ directory
│   └──...
├── mechanism
│   └──...
├── .editorconfig
├── CONTRIBUTING.md
├── README.md
├── SIMULATION_SETTINGS.md
└── try_excitations.ipynb
```

 > If you add a new control parameter, make sure it is in: documentation, Python and Matlab interfaces, example JSONs, ControlParameters and ParameterCombinator header and source files, including constructor, json parsing, string conversion.



## How it works?

This code is object-oriented. Most classes include some of the following methods. `ControlParameters` from [./include/control_parameters.h](./include/control_parameters.h) is an excellent example.
* **constructor with a builder**: many classes have a builder struct, like `ControlParameters::Builder`. This struct allows the usage of designated initializer lists. You no longer have to know the order of constructor arguments, and default values are handled as well.
* **constructor from JSON**: some classes can be constructed from a JSON to use with the interface. They accept both the file paths for the JSON and `nlohmann::ordered_json`. Some classes also have a `to_json()` method to help save as file.
* **to_string()**: use this to convert the object to a human readable string. Some arguments might be available, such as `colored` to print with some syntax coloring. Or `with_code`, which helps to print objects as valid C++ code (as builder structs). This is done in the form of builder struct, which you can copy into your own code.
* **ostream operator<< overload**: use to print objects directly through `std::cout`.
* **to_csv()**: convert the object to a CSV (string). Member variables are printed, separated by a comma. Order is defined by `csv_header`.

### High level overview

#### Running a single simulation:

To facilitate the calculations, two data container objects are required. The first is an instance of the `ControlParameters` class, which stores user-defined settings. These control parameters define the initial conditions, ambient properties, material properties, and acoustic excitation. They stay constant during the simulation. 
 
One of the control parameters is the name of the reaction mechanism to be used. This is translated into a pointer referencing an instance of the Parameters class. Multiple simulations with different control parameters can share the mechanism coefficients. These coefficients are not hard-coded at compile time but are imported at runtime from comment-augmented JSON files located at the [./mechanism](./mechanism) directory. The parsing overhead is negligible and occurs only once per mechanism, even in parameter studies. 

The  `ControlParameters` and `Parameters` classes serve exclusively as data containers and perform no computations. They are passed into an instance of the `OdeFun` class via the `init()` method, which also allocates memory for temporary results and precomputes certain values. The `OdeFun` object then sets up the ODE solver by calling the `initial_conditions()` method, which is executed only once. Subsequently, the right-hand-side function can be evaluated through the overload of `operator()`, the operator for calling the object as a function. This is the most computationally heavy portion of the simulation and is invoked at least ten thousand times by the ODE solver. 

The ODE solver is based on the SUNDIAL CVODE library and is implemented in the `OdeSolver` class. The `solve()` method takes an `OdeFun` instance as input and returns the numerical solution. The solution contains the time series, control parameters, solver metadata, and post-processing results. 


When using the interface, the control parameters (inputs) are read from a JSON file, and the numerical results (outputs) are written back into the same file. 

```
                                                 ┌─────────────────────────────────────┐
            ┌────────────────────────────────┐   │JSON Mechanism                       │
            │Parameters                      │   ├─────────────────────────────────────┤
            ├────────────────────────────────┤   │Constants and coefficients describing│
            │Holds constants and coefficients│   │reaction mechanisms, like:           │
            │describing reaction mechanisms. │---│∘ chemkin_elte2016_hydrogen          │
            │                                │   │∘ chemkin_kaust2023_ammonia          │
            │Parsed from JSON files.         │   │∘ chemkin_otomo2018_ammonia          │
            └────────────────────────────────┘   │∘ ...                                │
                            |                    └─────────────────────────────────────┘
                            |                                                          
                            |                                                          
              ┌──────────────────────────────┐                                           
┌──────────┐  │ControlParameters             │                                           
│JSON      │  ├──────────────────────────────┤                                           
├──────────┤--│Holds all settings influencing│                                           
│Input file│  │the simulation: R_E, P_amb,   │                                           
└──────────┘  │excitation parameters, ...    │                                           
              └──────────────────────────────┘                                           
                            |                                                          
┌─────────────────────────────────────────────────────┐                               
│OdeFun                                               │                               
├─────────────────────────────────────────────────────┤                               
│Computes the right-hand-side function: dxdt = f(x, t)│                               
│--                                                   │                               
│+ init(ControlParameters& cpar)                      │                               
│+ initial_conditions(double* x)                      │                               
│+ operator(double t, double* x, double* dxdt)        │                               
└─────────────────────────────────────────────────────┘                               
                            |                                                          
                                                                                        
            ┌──────────────────────────────────┐                                         
            │OdeSolver                         │                                         
            ├──────────────────────────────────┤                                         
            │Uses SUNDIALS CVODE to compute the│                                         
            │numerical solution.               │                                         
            │--                                │                                         
            │+ solve(OdeFun* ode_ptr, ...)     │                                         
            └──────────────────────────────────┘                                         
                            |                                                          
        ┌───────────────────────────────────────┐                                      
        │SimulationData                         │                                      
        ├───────────────────────────────────────┤   ┌───────────┐                      
        │Contains the results of the simulation:│   │JSON       │                      
        │∘ ControlParameters: inputs            │   ├───────────┤                      
        │∘ OdeSolution: numerical solution      │---│Output file│                      
        │∘ post-processing data: output         │   └───────────┘                      
        │--                                     │                                      
        │+ postprocess()                        │                                      
        └───────────────────────────────────────┘                                               
```

#### Running a parameter study

Investigating the dependence of simulation outputs on control parameters, or optimizing these outputs, requires executing the simulation with multiple sets of control parameters. The simplest approach is a brute‑force parameter study. In this method, certain parameters are held constant while others are subdivided with a specified resolution. The number of simulations grows exponentially with the number of parameters considered. For example, a 2D parameter study with a resolution of $100$ in each dimensions requires $10^4$ simulations, whereas a 4D study with the same resolution requires $10^8$.

Brute‑force parameter studies in the C++ implementation exploit multiple CPU threads to reduce runtime. Thread management is handled by the `ParameterStudy` class, with each thread executing an independent simulation in serial, similar to the workflow shown in [Running a single simulation](#running-a-single-simulation). In this case, the full time series is not retained, only the first and last time steps are saved. Each thread writes its results to a dedicated CSV file, with each row corresponding to one simulation. Upon completion, the `ParameterCombinator` class provides the next set of control parameters. A parameter study produces a directory containing CSV files, error logs, console output logs, and metadata (including hardware specifications, runtime, and success rate). For reproducibility, the study settings are saved both in JSON format (for use with the interface) and as builder structs (for use directly in code). 

```
                                    ┌──────────┐                                    
                                    │JSON      │                                    
                                    ├──────────┤                                    
                                    │Input file│                                    
                                    └──────────┘                                    
                                        |                                         
                                        |                                         
                    ┌─────────────────────────────────────┐                       
                    │ParameterCombinator                  │                       
                    ├─────────────────────────────────────┤                       
                    │Generates parameter combinations.    │                       
                    │--                                   │                       
                    │+ get_total_combination_count() const│                       
                    │+ get_next_combination()             │                       
                    └─────────────────────────────────────┘                       
                                        |                                         
    ┌────────────────────────────────────────────────────────────────────┐       
    │ParameterStudy                                                      │       
    ├────────────────────────────────────────────────────────────────────┤       
    │Runs parameter studies, distributing tasks across available threads.│       
    │--                                                                  │       
    │+ run()                                                             │       
    └────────────────────────────────────────────────────────────────────┘       
        |                    |                    |                    |            
        |                    |                    |                    |           
┌─────────────────┐  ┌─────────────────┐   ┌─────────────────┐   ┌─────────────────┐
│ControlParameters│  │ControlParameters│   │ControlParameters│   │ControlParameters│
├─────────────────┤  ├─────────────────┤   ├─────────────────┤   ├─────────────────┤
└─────────────────┘  └─────────────────┘   └─────────────────┘   └─────────────────┘
        |                     |                     |                    |         
    ┌──────┐              ┌──────┐              ┌──────┐             ┌──────┐      
    │OdeFun│              │OdeFun│              │OdeFun│             │OdeFun│      
    ├──────┤              ├──────┤              ├──────┤             ├──────┤      
    └──────┘              └──────┘              └──────┘             └──────┘      
        |                     |                     |                    |         
    ┌─────────┐          ┌─────────┐           ┌─────────┐           ┌─────────┐    
    │OdeSolver│          │OdeSolver│           │OdeSolver│           │OdeSolver│    
    ├─────────┤          ├─────────┤           ├─────────┤           ├─────────┤    
    └─────────┘          └─────────┘           └─────────┘           └─────────┘    
        |                     |                     |                    |         
        |                     |                     |                    |         
   ┌───────────┐        ┌───────────┐         ┌───────────┐         ┌───────────┐   
   │OdeSolution│        │OdeSolution│         │OdeSolution│         │OdeSolution│   
   ├───────────┤        ├───────────┤         ├───────────┤         ├───────────┤   
   └───────────┘        └───────────┘         └───────────┘         └───────────┘   
           |                     |                     |                    |         
   ┌────────────────┐    ┌────────────────┐    ┌────────────────┐   ┌────────────────┐ 
   │CSV             │    │CSV             │    │CSV             │   │CSV             │ 
   ├────────────────┤    ├────────────────┤    ├────────────────┤   ├────────────────┤ 
   │Thread 1 results│    │Thread 2 results│    │Thread 3 results│   │Thread 4 results│ 
   └────────────────┘    └────────────────┘    └────────────────┘   └────────────────┘ 
```

### Error handling

See `Error` and `ErrorHandler` classes in [./include/common.h](./include/common.h). Each error has:
* A severity (info, warning, error) and a type (general, preprocessing, ...). The type is meant to be used for statistical handling or simulation failiures.
* A string error message, usually with detailed context.
* Further metadata: exact timedate, location (file, function, line), and ID of the related control parameters combination, if the error is related to a specific simulation.

Errors are logged with the `LOG_ERROR(...)` macro, e.g.: `size_t error_ID = LOG_ERROR(Error::severity::error, Error::type::preprocess, message_string);`. Some classes have a member called `error_ID` used to store and propagate a state of error. `ErrorHandler::no_error` is used if the object is functioning as expected. Errors are printed as:

```
2025.02.15 21:31:08: ERROR (general) here is the error message with details in ./src/some_file.cpp:42: Class::function()
```

`ErrorHandler` is thread-safe. Use `ErrorHandler::set_log_file()` to log into a file, and set `ErrorHandler::print_when_log = false;` to disable printing on the console.




### Control parameters, and ODE function

The `ControlParameters` class, declared in [./include/control_parameters.h](./include/control_parameters.h), holds all parameters, which may influence the simulation. The recommended way to initialize is with builder struct and a designated initializer list: 
 
 ```cpp
 ControlParameters cpar = ControlParameters{ControlParameters::Builder{
     .ID                          = 0,
     .mechanism                   = "chemkin_elte2016_hydrogen",
     .R_E                         = 1.00000000000000008e-05,    // bubble equilibrium radius [m]
     .ratio                       = 1.00000000000000000e+00,    // R_0/R_E for unforced oscillations [-]
     .species                     = {"O2"},                     // names of species in initial bubble (array of strings)
     .fractions                   = {1.00000000000000000e+00},  // molar fractions of species in initial bubble (array of doubles)
     .P_amb                       = 1.01325000000000000e+05,    // ambient pressure [Pa]
     .T_inf                       = 2.93149999999999977e+02,    // ambient temperature [K]
     .alpha_M                     = 3.49999999999999978e-01,    // water accommodation coefficient [-]
     .P_v                         = 2.33809999999999991e+03,    // vapour pressure [Pa]
     .mu_L                        = 1.00000000000000002e-03,    // dynamic viscosity [Pa*s]
     .rho_L                       = 9.98200000000000045e+02,    // liquid density [kg/m^3]
     .c_L                         = 1.48300000000000000e+03,    // sound speed [m/s]
     .surfactant                  = 1.00000000000000000e+00,    // surface tension modifier [-]
     .enable_heat_transfer        = true,
     .enable_evaporation          = true,
     .enable_reactions            = true,
     .enable_dissipated_energy    = true,
     .enable_van_der_waals        = true,
     .enable_gilmore              = true,
     .enable_nasg                 = true,
     .enable_rate_thresholding    = true,
     .target_specie               = "H2",
     .excitation_type             = Parameters::excitation::sinusoid,                      // type of excitation (enum)
     .excitation_params           = {-2.00000000000000000e+05, 3.00000000000000000e+04},   // parameters for excitation (array of doubles)
     .excitation_cycles           = 1,                                                     // number of excitation cycles to use (according to freq/freq1 in excitation_params) [-]
     .ramp_up_cycles              = 0                                                      // number of cycles until the excitation reaches full amplitude (0<=ramp_up_cycles<=excitation_cycles/2) [-]
 }};
 ```
 These are also the default values, any builder argument can be missed. Use the `to_string()` method or the `ostream operator<<` overload to print the control parameters to the console. It will be printed in the same format as above, which is also valid code. The class also has a CSV header and a `to_csv()` method, like many other classes.
 
 The right-hand side of the ODE is calculated by the `OdeFun` class declared in [./include/ode_fun.h](./include/ode_fun.h). (The right-hand side is $f()$ if the ODE is $\frac{dx}{dt} = f(x)$ and the initial condition is $x(0) = x_0$.) `OdeFun` is a function object, with the call operator (`operator()`) defined. It has to be initialized, so it can allocate and reuse arrays for temporary results, and not call `new` during high performance calculations:
 
 ```cpp
 OdeFun ode;
 ode.init(cpar);
 ```
 
 `OdeFun` is also responsible to calculate the initial conditions, however, an array of appropriate size has to be provided as a pointer. This step is handled by the solver:
 
 ```cpp
 std::vector<double> x_0(cpar.par->num_species+4);
 ode.initial_conditions(x_0.data());
 ```
 
 Now, the right-hand side can be calculated with `operator()`. The pointers for $x$ and $\frac{dx}{dt}$ have to be provided. This function returns the type `is_success`, which is effectively a boolean. However, in case of an error, `ode.cpar.error_ID` will also change from `Error::no_error`. You may obtain the error as:
 ```cpp
 Error error = ErrorHandler.get_error(ode.cpar.error_ID);
 std::cout << error <<std::endl;
 ```
 
 ### Running a simulation
 
 ```cpp
 // init cpar and ode
 ControlParameters cpar = ControlParameters{{ /* ... */ }};
 OdeFun ode;
 ode.init(cpar);
 
 // solve the ODE      
 OdeSolver solver(cpar.par->num_species+4);
 SimulationData data = solver.solve(
     1.0,     // t_max [s]
     &ode,    // ode_ptr
     60.0,    // timeout [s]
    true     // whether to save solution (or just first and last step)
 );
 
 std::cout << data << std::endl;
 ```
 
 ### Running a bruteforce parameter study
 
 A brute-force parameter study can be defined by the `ParameterCombinator` class declared in [./include/parameter_combinator.h](./include/parameter_combinator.h). It is also initialized with a builder struct and a designated initializer list, similar to `ControlParameters`. However, arguments have to be children of the `Range` class: `Const`, `LinearRange`, `LogRange`, `GeomRange`. It also has defaults for missing arguments, and can be printed to console in usable code format:
 
 ```cpp
 ParameterCombinator parameter_combinator = ParameterCombinator{ParameterCombinator::Builder{
     .mechanism                   = Parameters::mechanism::chemkin_elte2016_hydrogen,
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
     .excitation_type             = Parameters::excitation::sinusoid
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


 ### Testing and benchmarking with built in unit testing framework

A lightweight custom framework was developed for this project, see [./include/test.h](./include/test.h). This is a minimalist framework, constituted by only 160 lines of code. Tests are organized into groups and groups into test cases. Test cases may contain an arbitrary number of assertions (comparisons), defined by macros. For example:
```C++
testing::Tester tester{"Title of test group"}; 
 
ADD_TEST(tester, "Successful test case", 
    ASSERT_TRUE(1 + 1 == 2); 
); 
ADD_TEST(tester, "Failing test case", 
    ASSERT_TRUE(1 + 1 == 3); 
); 
 
tester.run_tests(); 
testing::Tester::print_summary(); 
```
```
Title of test group 
    PASSED: Successful test case 
    FAILED: Failing test case  
            Assert true failed: 1 + 1 == 3 != true in ./src/main.cpp:31 
_____________________________________________ 
TOTAL PASSED: 1
TOTAL FAILED: 1
```

Available assertions in the framework:

| Macro | Description |
|:---|:---|
| `ASSERT_TRUE(condition)` | Asserts that `condition` evaluates to `true`. |
| `ASSERT_FALSE(condition)` | Asserts that `condition` evaluates to `false`. |
| `ASSERT_EQUAL(a, b)` | Asserts that `a` is strictly equal to `b`. |
| `ASSERT_NOT_EQUAL(a, b)` | Asserts that `a` is not equal to `b`. |
| `ASSERT_NEAR(a, b, tol)` | Asserts that the absolute difference between `a` and `b` is $\le$ `tol`. |
| `ASSERT_APPROX(a, b, tol)` | Asserts that the relative difference between `a` and `b` is $\le$ `tol`. |
| `ASSERT_APPROX_ARRAY(a, b, size, tol)` | Checks `ASSERT_APPROX` for each element up to index `size`. |
| `FAIL(message)` | Unconditionally fails the test case with the specified `message`. |

This framework is well suited for floating point numbers and arrays:
```C++
ADD_TEST(tester, "Compare arrays", 
    std::vector<double> reference_data = {1.0, 2.0, 3.0}; 
    auto test_data = function_to_test();    // returns {1.0, 2.0, 3.1} 
    ASSERT_APPROX_ARRAY(reference_data, test_data, test_data.size(), 1e-10); 
); 
```
```
FAILED: Compare arrays 
        Assert approx failed: reference_data[2] != test_data[2] (3.000000 != 3.100000) in ./src/main.cpp:35 
```

For performance testing, the project provides a simple benchmarking macro in [./include/common.h](./include/common.h). You provide the line of code to test, and the number of repetitions `N`. It will run the statement `N` times and print the average runtime in microseconds/milliseconds to standard output: 

```C++
BENCHMARK_LINE(function_to_benchmark(1.14), 1000);
```
```
Runtime of function_to_benchmark(1.14) is 35.72 us 
```

**Where to define tests?**
All tests and benchmarks should be defined in source files located under the [./test/](./test/) directory. To run tests, compile with the `--test` flag. To run benchmarks, compile with the `--benchmark` flag (e.g., `python build.py --run --benchmark --optimize2`).
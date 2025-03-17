# Bubble dynamics simulation C++
A C++ powered version of [Bubble_dynamics_simulation](https://github.com/kozakaron/Bubble_dynamics_simulation).

## Prerequisits

 * Python3 (with numpy, pygments)
 * gnu g++/llvm clang++ compiler
 * gnu gdb/llvm lldb debugger

 ## Compilation

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
 The script collects each *.cpp* file automatically, and compiles them into object files, each one on a seperate thread. Then the object files are linked togeather. Syntax colored error messages are printed grouped by *.cpp* file. Use flags `--test`, `--benchmark` to run unit tests or benchmarks respectively. Use `--run` to run the executable automatically after compilation. 

For example you may compile+run your tests and benchmarks as: `python ./dev/build.py --run --warning --optimize3 --benchmark --test` or `python ./dev/build.py -r -w -o3 -b -t`.

Notes:
 * You may have problem with the **sanitizers**, especially on Windows. Feel free to turn them off. Also turn them off for high performance computing. However, they are very powerful tools to find bugs. I never managed to make the thread sanitizer work.
 * Some systems strugle with the long double. If the long double exists unit test fails or you encounter errors, you may try to change `-mlong-double-80` to 64 or 128. (I actually found a compiler error, causing segfault while calculating with long doubles on windows: [link](https://github.com/llvm/llvm-project/issues/118138))
 * Thanks to multithreading, compilation performance is very good. However, in the future, it can be further improved. First, by using precompiled headers. And second, by creating a hash of each preprocessed source file, and only recompileing the ones that changed.

## How it works?
### Error handling

See `Error` and `ErrorHandler` classes in [./include/common.h](./include/common.h). Each error has:
 * A severity (info, warning, error) and a type (general, preprocessing, ...). The type is meant to be used for statistical handling or simulation failiures.
 * A string error meassage, usually with detailed context.
 * Further metadata: exact timedate, location (file, function, line), and ID of the related control parameters combination, if the error is related to a simulation.

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

`OdeFun` is also responsible to calculate the initial conditions, however, an array of appropiate size has to be provided as a pointer:

```cpp
std::vector<double> x_0(ode.par->num_species+4);
ode.initial_conditions(x_0.data());
```

Now, the right hand side can be calculated with `operator()`. The pointers for $x$ and $\frac{dx}{dt}$ has to be provided. ODE solvers usually take a function, thus, a lambda has to be defined:

```cpp
auto ode_fun = [&ode](const double t, const double *x, double *dxdt) -> is_success { return ode(t, x, dxdt); };
```

The type `is_success` is a boolean. However, in case of an error, `ode.cpar.error_ID` will also change from `Error::no_error`. You may obtain the error as:
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

// calculate initial conditions, turn ode function object into ode_fun lambda function.
std::vector<double> x_0(ode.par->num_species+4);
ode.initial_conditions(x_0.data());
auto ode_fun = [&ode](const double t, const double *x, double *dxdt) -> is_success { return ode(t, x, dxdt); };

// solve the ODE
RKCK45 solver;
solver.solve(0.0, t_max, (double*)x_0.data(), ode.par->num_species+4, ode_fun, &ode.cpar.error_ID, 60.0, false);

// get solution, do postprocessing
OdeSolution sol = solver.get_solution();
SimulationData data(cpar, sol);
std::cout << data << std::endl;
```

### Running a bruteforce parameter study

A bruteforce parameter study can be defined by the `ParameterCombinator` class declared in [./include/parameter_study.h](./include/parameter_study.h). It is also initialized with a builder struct and a designated initializer list, similar to `ControlParameters`. However, arguments have to be a childre of the `Range` class: `Const`, `LinearRange`, `PowRange`. It also has defaults for missing arguments, and can be printed to console in usable code format:

```cpp
ParameterCombinator parameter_study = ParameterCombinator{ParameterCombinator::Builder{
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

std::cout << parameter_study << std::endl;
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
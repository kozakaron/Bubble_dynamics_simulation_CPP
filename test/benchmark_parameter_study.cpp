#ifdef BENCHMARK

#include "common.h"
#include "parameters.h"
#include "control_parameters.h"
#include "ode_fun.h"
#include "ode_solver.h"
#include "parameter_study.h"
#include "test_list.h"

#include <thread>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <mutex>

std::mutex mutex;
const double t_max = 1000.0e-6;
const size_t num_threads = std::thread::hardware_concurrency() - 1;
double best_energy_demand = SimulationData::infinite_energy_demand;
ParameterCombinator parameter_combinator = ParameterCombinator{ParameterCombinator::Builder{
    .mechanism                   = Parameters::mechanism::chemkin_ar_he,
    .R_E                         = LinearRange(0.000005, 0.000125, 50),                // {5e-06, 7.44898e-06, 9.89796e-06, 1.23469e-05, 1.47959e-05, 1.72449e-05, 1.96939e-05, ..., 0.000122551, 0.000125}
    .species                     = {"AR"},
    .fractions                   = {1.00000000000000000e+00},
    .P_amb                       = Const(101325.000000),                               // {101325}
    .T_inf                       = Const(293.150000),                                  // {293.15}
    .alfa_M                      = Const(0.350000),                                    // {0.35}
    .P_v                         = Const(2338.100000),                                 // {2338.1}
    .mu_L                        = Const(0.001000),                                    // {0.001}
    .rho_L                       = Const(998.200000),                                  // {998.2}
    .c_L                         = Const(1483.000000),                                 // {1483}
    .surfactant                  = Const(1.000000),                                    // {1}
    .enable_heat_transfer        = true,
    .enable_evaporation          = true,
    .enable_reactions            = true,
    .enable_dissipated_energy    = true,
    .target_specie               = "H2",
    .excitation_params           = {
        LinearRange(-100000.000000, -300000.000000, 50),    // {-100000, -104082, -108163, -112245, -116327, -120408, -124490, ..., -295918, -300000}
        Const(20000.000000),                                // {20000}
        Const(1.000000)                                     // {1}
    },
    .excitation_type             = Parameters::excitation::sin_impulse
}};


void parameter_study_task()
{
    OdeFun ode;
    RKCK45 solver;
    auto ode_fun = [&ode](const double t, const double *x, double *dxdt) -> is_success { return ode(t, x, dxdt); }; 

    while(true)
    {
        auto [success, cpar] = parameter_combinator.get_next_combination();
        if (!success) break;

        ode.init(cpar);
        std::vector<double> x_0(ode.par->num_species+4);
        ode.initial_conditions(x_0.data());

        solver.solve(0.0, t_max, (double*)x_0.data(), ode.par->num_species+4, ode_fun, &ode.cpar.error_ID, 60.0, false);
        OdeSolution sol = solver.get_solution();
        SimulationData data(cpar, sol);

        {
            std::unique_lock<std::mutex> lock(mutex);
            best_energy_demand = std::min(best_energy_demand, data.energy_demand);
            std::cout << data.to_small_string(parameter_combinator, best_energy_demand, true) << "          \r";
            // \r will allow the next line to overwrite this one (if the line fits in the terminal) 
        }
    }
}

void benchmark_parameter_study()
{
    std::cout << colors::bold << "Small parameter study with RKCK45 solver and chemkin_ar_he mechanism" << colors::reset << std::endl;
    std::cout << "total_combination_count = " << parameter_combinator.get_total_combination_count() << std::endl;
    //std::cout << parameter_combinator << std::endl;

    Timer timer;
    timer.start();
    std::vector<std::thread> threads(num_threads);

    for (size_t i = 0; i < num_threads; i++)
    {
        threads[i] = std::thread(parameter_study_task);
    }

    for (size_t i = 0; i < num_threads; i++)
    {
        threads[i].join();
    }
    

    double runtime = timer.lap();
    std::cout << "\n\nTotal runtime: " << Timer::format_time(runtime) << std::endl;
    std::cout << "\nAverage runtime per combination: " << Timer::format_time(runtime / parameter_combinator.get_total_combination_count()) << std::endl;
}

#endif // BENCHMARK
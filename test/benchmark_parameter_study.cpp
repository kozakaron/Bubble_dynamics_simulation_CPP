#ifdef BENCHMARK
#include <thread>
#include <iostream>
#include <string>

#include "parameter_study.h"
#include "ode_solver_sundials.h"


std::string save_folder_base_name = "./_parameter_studies/test";
const double t_max = 1.0;
const double timeout = 60.0;
const size_t num_threads = std::thread::hardware_concurrency();

ParameterCombinator parameter_combinator = ParameterCombinator{ParameterCombinator::Builder{
    .mechanism                   = Parameters::mechanism::chemkin_ar_he,
    .R_E                         = LinearRange(0.000005, 0.000125, 20),                // {5e-06, 1.13158e-05, 1.76316e-05, 2.39474e-05, 3.02632e-05, 3.65789e-05, 4.28947e-05, ..., 0.000118684, 0.000125}
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
        LinearRange(-100000.000000, -300000.000000, 20),    // {-100000, -110526, -121053, -131579, -142105, -152632, -163158, ..., -289474, -300000}
        Const(20000.000000),                                // {20000}
        Const(1.000000)                                     // {1}
    },
    .excitation_type             = Parameters::excitation::sin_impulse
}};


void benchmark_parameter_study()
{
    std::cout << colors::bold << "Small parameter study with SUNDIALS CVODE solver and chemkin_ar_he mechanism" << colors::reset << std::endl;
    std::cout << "total_combination_count = " << parameter_combinator.get_total_combination_count() << std::endl;
    
    ErrorHandler::print_when_log = true;
    ParameterStudy parameter_study = ParameterStudy{
        parameter_combinator,
        save_folder_base_name,
        t_max,
        timeout
    };
    parameter_study.run(num_threads, true);
    ErrorHandler::clear_errors();
}

#endif // BENCHMARK
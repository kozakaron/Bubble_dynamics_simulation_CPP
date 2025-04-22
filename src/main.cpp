#include <thread>

#include "common.h"
#include "parameters.h"
#include "control_parameters.h"
#include "ode_fun.h"
#include "ode_solver.h"
#include "ode_solver_sundials.h"
#include "parameter_study.h"
#include "test_list.h"

#include "nlohmann/json.hpp"
#include "cxxopts.hpp"


using namespace std;
using namespace nlohmann;


OdeSolver* solver_factory(size_t num_dim)
{
    return new OdeSolverCVODE(num_dim);
}


int main(int argc, char **argv)
{
    // Initialize logging
    cxxopts::ParseResult result;
    cxxopts::Options options("Bubble dynamics simulation C++", "A high performance SUNDIALS CVODE powered zero dimensional (ODE) sonochemistry simulation of an acustically excited bubble.");
    options.add_options()
        ("h,help", "Show help")
        ("v,version", "Show version")
        ("r,run", "Run a simulation with the given JSON file", cxxopts::value<std::string>())
        ("tmax", "Simulation end time in seconds", cxxopts::value<double>()->default_value("1.0"))
        ("timeout", "Timeout in seconds", cxxopts::value<double>()->default_value("60.0"))
        ("save", "Set this flag to save all timesteps, skip it to save only the first and last steps", cxxopts::value<bool>()->default_value("false"))
        ("log", "Set log file", cxxopts::value<std::string>())
        ("parameter_study", "Run a parameter study with the given JSON file", cxxopts::value<std::string>())
        ("directory", "Set save directory for parameter_study", cxxopts::value<std::string>()->default_value("./_parameter_studies/test"))
        ;
    
    // Parse command line arguments
    try
    {
        result = options.parse(argc, argv);
    }
    catch(const std::exception& e)
    {
        LOG_ERROR("Error parsing command line arguments: " + std::string(e.what()));
        return 1;
    }
    
    // Print help
    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        return 0;
    }

    // Print version
    if (result.count("version"))
    {
        std::cout << "Version: " << VERSION << std::endl;
        return 0;
    }

    // Init logging
    if (result.count("log"))
    {
        std::string log_path = result["log"].as<std::string>();
        ErrorHandler::set_log_file(log_path);
    }

    // Run simulation
    if (result.count("run"))
    {
        std::string json_path = result["run"].as<std::string>();

        // Run simulation
        ControlParameters cpar(json_path);        
        OdeFun ode; ode.init(cpar);
        OdeSolverCVODE solver(ode.par->num_species+4);
        OdeSolution solution = solver.solve(
            result["tmax"].as<double>(),    // t_max [s]
            &ode,                           // ode_ptr
            result["timeout"].as<double>(), // timeout [s]
            result["save"].as<bool>()       // save solution
        );
        SimulationData data(cpar, solution);

        // Save results
        data.save_json_with_binary(json_path);
    }

    // Run parameter study
    else if (result.count("parameter_study"))
    {
        std::string json_path = result["parameter_study"].as<std::string>();
        ParameterCombinator parameter_combinator(json_path);
        if (ErrorHandler::get_error_count() != 0) return 1;
        const size_t num_threads = std::thread::hardware_concurrency();
        ParameterStudy parameter_study(
            parameter_combinator,
            result["directory"].as<std::string>(),  // save_folder_base_name
            solver_factory,                         // solver_factory
            result["tmax"].as<double>(),            // t_max [s]
            result["timeout"].as<double>()          // timeout [s]
        );
        parameter_study.run(num_threads, true);
    }
    
    
#ifdef TEST
    testing::test_common();
    testing::test_par_cpar();
    testing::test_ode_fun_otomo2018();
    testing::test_ode_fun_ar_he();
    testing::test_ode_solver();
    testing::test_parameter_study();
    testing::print_test_summary();
#endif  // TEST

#ifdef BENCHMARK
    benchmark_ode_fun();
    benchmark_speedup();
    benchmark_parameter_study();
#endif  // BENCHMARK

    if (ErrorHandler::get_error_count() != 0)  return 1;
    return 0;
}
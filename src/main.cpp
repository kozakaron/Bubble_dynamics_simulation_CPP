#include "common.h"
#include "parameters.h"
#include "control_parameters.h"
#include "ode_fun.h"
#include "ode_solver.h"
#include "ode_solver.h"
#include "parameter_study.h"
#include "test_list.h"

#include "nlohmann/json.hpp"
#include "cxxopts.hpp"

using namespace std;
using namespace nlohmann;


int main(int argc, char **argv)
{
    // ControlParameters cpar=ControlParameters(ControlParameters::Builder{
		// .ID                          = 0,
		// .mechanism                   = Parameters::mechanism::chemkin_kaust2023_ammonia_oxygenless,
		// .R_E                         = 1.00000000000000008e-05,    // bubble equilibrium radius [m]
		// .species                     = {"O2","AR"},
		// .fractions                   = {0.73,0.27},
		// .P_amb                       = 1.01325000000000000e+05,    // ambient pressure [Pa]
		// .T_inf                       = 2.93149999999999977e+02,    // ambient temperature [K]
		// .alfa_M                      = 3.49999999999999978e-01,    // water accommodation coefficient [-]
		// .P_v                         = 2.33809999999999991e+03,    // vapour pressure [Pa]
		// .mu_L                        = 1.00000000000000002e-03,    // dynamic viscosity [Pa*s]
		// .rho_L                       = 9.98200000000000045e+02,    // liquid density [kg/m^3]
		// .c_L                         = 1.48300000000000000e+03,    // sound speed [m/s]
		// .surfactant                  = 1.00000000000000000e+00,    // surface tension modifier [-]
		// .enable_heat_transfer        = true,
		// .enable_evaporation          = true,
		// .enable_reactions            = true,
		// .enable_dissipated_energy    = true,
		// .target_specie               = "H2",
		// .excitation_params           = {-2.00000000000000000e+05, 3.00000000000000000e+04, 1.00000000000000000e+00},
		// .excitation_type             = Parameters::excitation::sin_impulse}
	// );
	// cpar.R_E =1e-4;
	// cout << cpar << endl;
	
	// OdeFun ode;
	// ode.init(cpar);
	// cout << ode.par->model <<endl;
	
	// OdeSolverCVODE solver(ode.par->num_species+4);
	// OdeSolution sol= solver.solve(1.0,&ode);
	// cout << sol <<endl;
	
	// SimulationData data = SimulationData(cpar,sol);
	// cout << data <<endl;
	
	// ErrorHandler::set_log_file("error.log");
	
	
	// size_t errorid = LOG_ERROR(Error::severity::info,Error::type::general,"Kedves LOG_ERROR! Szia!");
	// Error error = ErrorHandler::get_error(errorid);
	// cout << error <<endl;
	// (void)ErrorHandler::no_error;
	
	
	// return 0;
	
	
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
        ("save_jacobian", "Set this flag to save the evaluated Jacobian matrixes", cxxopts::value<bool>()->default_value("false"))
        ("log", "Set log file", cxxopts::value<std::string>())
        ("parameter_study", "Run a parameter study with the given JSON file", cxxopts::value<std::string>())
        ("directory", "Set save directory for parameter_study", cxxopts::value<std::string>()->default_value("./_parameter_studies/test"))
        ("cpu", "Set number of CPU cores to use for parameter_study", cxxopts::value<size_t>()->default_value(std::to_string(std::thread::hardware_concurrency())))
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
        OdeSolver solver(cpar.par->num_species+4);
        SimulationData data = solver.solve(
            result["tmax"].as<double>(),        // t_max [s]
            &ode,                               // ode_ptr
            result["timeout"].as<double>(),     // timeout [s]
            result["save"].as<bool>(),          // save solution
            result["save_jacobian"].as<bool>()  // save jacobian
        );

        // Save results
        data.save_json_with_binary(json_path);
    }

    // Run parameter study
    else if (result.count("parameter_study"))
    {
        std::string json_path = result["parameter_study"].as<std::string>();
        ParameterCombinator parameter_combinator(json_path);
        if (ErrorHandler::get_error_count() != 0) return 1;
        size_t num_threads = result["cpu"].as<size_t>();
        ParameterStudy parameter_study(
            parameter_combinator,
            result["directory"].as<std::string>(),  // save_folder_base_name
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
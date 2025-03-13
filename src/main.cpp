#include "common.h"
#include "parameters.h"
#include "control_parameters.h"
#include "ode_fun.h"
#include "ode_solver.h"
#include "ode_solver_sundials.h"
#include "parameter_study.h"
#include "test_list.h"


using namespace std;

int main(int argc, char **argv)
{
    (void)argc; (void)argv;
    //ErrorHandler::set_log_file("log.txt");

// solve with cvode
    ControlParameters cpar = ControlParameters{{ 
        .mechanism = Parameters::mechanism::chemkin_otomo2018,
        .species = {"H2", "N2"},
        .fractions = {0.75, 0.25},
    }};
    OdeFun ode; ode.init(cpar);
    OdeSolverCVODE cvode(ode.par->num_species+4);
    cout << cpar << endl;

    OdeSolution solution = cvode.solve(1.0, &ode, 2, false);

    cout << solution << endl;
    
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

    //ErrorHandler::print_errors();
    if (ErrorHandler::get_error_count() != 0)  return 1;
    return 0;
}
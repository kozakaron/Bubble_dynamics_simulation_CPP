#include "common.h"
#include "parameters.h"
#include "ode_fun.h"
#include "test_list.h"


using namespace std;


int main(int argc, char **argv)
{
    (void)argc; (void)argv;
    //ErrorHandler::set_log_file("log.txt");

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
#endif  // BENCHMARK

    ErrorHandler::print_errors();
    if (ErrorHandler::get_error_count() != 0)  return 1;
    return 0;
}
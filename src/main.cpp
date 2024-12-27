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
    test_common();
    test_par_cpar();
    test_ode_fun();
    print_test_summary();
#endif  // TEST

#ifdef BENCHMARK
    benchmark_ode_fun();
#endif  // BENCHMARK

    ErrorHandler::print_errors();
    if (ErrorHandler::get_error_count() != 0)  return 1;
    return 0;
}
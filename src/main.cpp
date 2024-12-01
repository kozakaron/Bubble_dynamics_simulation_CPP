#include "common.h"
#include "parameters.h"
#include "ode_fun.h"
#include "test_list.h"


using namespace std;


int main(int argc, char **argv)
{
    (void)argc; (void)argv;

#ifdef TEST
    test_ode_fun();
#endif  // TEST

#ifdef BENCHMARK
    benchmark_ode_fun();
#endif  // BENCHMARK

    return 0;
}
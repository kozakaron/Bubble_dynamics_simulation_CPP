#ifndef TEST_LIST_H
#define TEST_LIST_H

#include "test.h"


#ifdef TEST
#ifndef CHEMKIN_OTOMO2018
#error CHEMKIN_OTOMO2018 must be defined for testing
#endif
void test_ode_fun();



#endif  // TEST


#ifdef BENCHMARK
#ifndef CHEMKIN_OTOMO2018
#error CHEMKIN_OTOMO2018 must be defined for benchmarking
#endif
void benchmark_ode_fun();



#endif  // BENCHMARK


#endif  // TEST_LIST_H
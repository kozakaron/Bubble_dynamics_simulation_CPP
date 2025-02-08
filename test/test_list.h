#ifndef TEST_LIST_H
#define TEST_LIST_H


#ifdef TEST

namespace testing{

void test_common();
void test_par_cpar();
void test_ode_fun_otomo2018();
void test_ode_fun_ar_he();
void test_ode_solver();
void test_parameter_study();


void print_test_summary();
    
}   // namespace testing
#endif  // TEST


#ifdef BENCHMARK
void benchmark_ode_fun();
void benchmark_speedup();



#endif  // BENCHMARK


#endif  // TEST_LIST_H
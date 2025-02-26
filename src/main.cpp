#include "common.h"
#include "parameters.h"
#include "control_parameters.h"
#include "ode_fun.h"
#include "ode_solver.h"
#include "parameter_study.h"
#include "test_list.h"

// libode
// git repo: https://github.com/markmbaum/libode
// in "libode\src\ode_adaptive.cc" change tolerances:
    // abstol_ = 1e-10;
    // reltol_ = 1e-10;
// in "libode\src\ode_io.cc" modify ode_print_exit(), and other places where exit(EXIT_FAILURE) is used:
    // throw std::runtime_error("ODE FAILURE: " + std::string(msg));
    // instead of exit(EXIT_FAILURE);
#include "ode_rkf_32.h"
#include "ode_rk_43.h"
#include "ode_rkck.h"
#include "ode_dopri_54.h"
#include "ode_vern_65.h"
#include "ode_vern_76.h"
#include "ode_dopri_87.h"
#include "ode_vern_98.h"

#include "ode_grk4a.h"
#include "ode_gauss_6.h"
#include "ode_lobatto_iiic_6.h"
#include "ode_radau_iia_5.h"
#include "ode_sdirk_43.h"

using namespace std;

#define PRINT_SOL(solvername, implicit)    cout << setw(strw) << solvername << setw(numw) << boolalpha << implicit << setw(numw) << Timer::format_time(timer.lap()) << setw(numw) << sys.get_nstep() << setw(numw) << sys.get_neval() << setw(numw) << sys.get_sol(ode.par->num_species+4) << endl;
#define SOLVE(solver) ode::OdeFun_libode<ode::solver> sys; \
    sys.ode = &ode; \
    sys.set_sol(0, 0.0); \
    for(size_t i = 0; i < ode.par->num_species+4; i++) \
        sys.set_sol(i+1, x_0[i]); \
    Timer timer; timer.start(); \
    sys.solve_adaptive(t_max, 1e-15, false); \

namespace ode {
template<class Integrator>
class OdeFun_libode : public Integrator {

    public:
        OdeFun *ode;

        OdeFun_libode () : Integrator (17) {}

        void ode_fun (double *solin, double *fout) {
            double t = solin[0];
            double* x = solin + 1;
            double* dxdt = fout + 1;
            fout[0] = 1.0;
            ode->operator()(t, x, dxdt);
        }
};
} // namespace ode

int main(int argc, char **argv)
{
    (void)argc; (void)argv;
    //ErrorHandler::set_log_file("log.txt");

// setup
    ControlParameters cpar;
    OdeFun ode;
    ode.init(cpar);
    std::vector<double> x_0(ode.par->num_species+4);
    ode.initial_conditions(x_0.data());
    double t_max = 1.0e-3;
    const int strw = 25;
    const int numw = 20;
    cout << "Test cpar: " << cpar << endl;
    cout << "Integration time interval end: " << Timer::format_time(t_max) << endl << endl;
    cout << setw(strw) << "solver name" << setw(numw) << "implicit" << setw(numw) << "runtime" << setw(numw) << "num_steps" << setw(numw) << "num_fun_eval" << setw(numw) << "x_last[-1]" << endl;

// solve with builtin RKCK45
    auto ode_fun = [&ode](const double t, const double *x, double *dxdt) -> is_success { return ode(t, x, dxdt); };
    RKCK45 solver;
    solver.solve(0.0, t_max, (double*)x_0.data(), ode.par->num_species+4, ode_fun, &ode.cpar.error_ID, 60.0, false);
    OdeSolution sol = solver.get_solution();
    cout << setw(strw) << "own RKCK45" << setw(numw) << boolalpha << false << setw(numw) << Timer::format_time(sol.runtime) << setw(numw) << sol.num_steps << setw(numw) << setw(numw) << sol.num_fun_evals << setw(numw) << sol.x.back().back() << endl;

// solve with libode
    ErrorHandler::print_when_log = false;
    {
        SOLVE(OdeRKF32)
        PRINT_SOL("libode RKF32", false);
    }
    {
        SOLVE(OdeRK43)
        PRINT_SOL("libode RK43", false);
    }
    {
        SOLVE(OdeRKCK)
        PRINT_SOL("libode OdeRKCK", false);
    }
    {
        SOLVE(OdeDoPri54)
        PRINT_SOL("libode OdeDoPri54", false);
    }
    try {
        SOLVE(OdeVern65)
        PRINT_SOL("libode OdeVern65", false);
    } catch (std::exception &e) {
        cout << setw(strw) << "libode OdeVern65" << setw(numw) << false << setw(numw) << "failed: " << e.what() << endl;
    }
    try {
        SOLVE(OdeVern76)
        PRINT_SOL("libode OdeVern76", false);
    } catch (std::exception &e) {
        cout << setw(strw) << "libode OdeVern76" << setw(numw) << false << setw(numw) << "failed: " << e.what() << endl;
    }
    try {
        SOLVE(OdeDoPri87)
        PRINT_SOL("libode OdeDoPri87", false);
    } catch (std::exception &e) {
        cout << setw(strw) << "libode OdeDoPri87" << setw(numw) << false << setw(numw) << "failed: " << e.what() << endl;
    }
    try {
        SOLVE(OdeVern98)
        PRINT_SOL("libode OdeVern98", false);
    } catch (std::exception &e) {
        cout << setw(strw) << "libode OdeVern98" << setw(numw) << false << setw(numw) << "failed: " << e.what() << endl;
    }
    try {
        SOLVE(OdeGRK4A)
        PRINT_SOL("libode GRK4A", true);
    } catch (std::exception &e) {
        cout << setw(strw) << "libode GRK4A" << setw(numw) << true << setw(numw) << "failed: " << e.what() << endl;
    }
    /*{
        SOLVE(OdeSDIRK43)
        PRINT_SOL("libode OdeSDIRK43", true);
    }*/
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

    ErrorHandler::print_errors();
    if (ErrorHandler::get_error_count() != 0)  return 1;
    return 0;
}
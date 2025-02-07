#ifdef TEST
#include <array>

#include "common.h"
#include "test_list.h"
#include "test.h"
#include "ode_solver.h"


NO_OPTIMIZATION

namespace testing{

using std::array;

size_t error_ID = ErrorHandler::no_error;
is_success ode_fun_hires(const double t, const double *x, double *dxdt)
{
    (void)t;
    dxdt[0] = -1.71*x[0] + 0.43*x[1] + 8.32*x[2] + 0.0007;
    dxdt[1] = 1.71*x[0] - 8.75*x[1];
    dxdt[2] = -10.03*x[2] + 0.43*x[3] + 0.035;
    dxdt[3] = 8.32*x[1] + 1.71*x[2] - 1.12*x[3];
    dxdt[4] = -1.745*x[4] + 0.43*x[5] + 0.43*x[6];
    dxdt[5] = -280*x[5]*x[7] + 0.69*x[3] + 1.71*x[4] - 0.43*x[5] + 0.69*x[6];
    dxdt[6] = 280*x[5]*x[7] - 1.81*x[6];
    dxdt[7] = -280*x[5]*x[7] + 1.81*x[6];
    return true;
}

void test_RKCK45();

void test_ode_solver()
{
    test_RKCK45();
}

void test_RKCK45()
{
    Tester tester = Tester("Test ode_solver.h's RKCK45 class");

    ADD_TEST(tester, "HIRES problem",
        array<const double, 8> x_0 = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0057};
        array<const double, 8> x_sol = {2.96363741e-02, 5.79179425e-03, 5.70765963e-03, 5.17391304e-02, 8.30071828e-01, 3.36411375e+00, 5.68906817e-03, 1.09318263e-05};
        RKCK45 solver;
        solver.solve(0.0, 321.8122, x_0.data(), 8, ode_fun_hires, &error_ID, 2.0);
        auto sol = solver.get_solution();
        ASSERT_TRUE(sol.success());
        ASSERT_NEAR(sol.t.back(), 321.8122, 1.0e-10);
        ASSERT_APPROX_ARRAY(sol.x.back().data(), x_sol.data(), 8, 1.0e-5);
    );

    ADD_TEST(tester, "HIRES problem with error",
        array<const double, 8> x_0 = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0057};
        RKCK45 solver;
        error_ID = 69;
        solver.solve(0.0, 321.8122, x_0.data(), 8, ode_fun_hires, &error_ID, 2.0);
        auto sol = solver.get_solution();
        ASSERT_FALSE(sol.success());
    );

    tester.run_tests();
}

}   // namespace testing

#endif // TEST
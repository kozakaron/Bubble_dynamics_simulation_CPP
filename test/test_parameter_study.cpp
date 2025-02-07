#ifdef TEST
#include <array>

#include "common.h"
#include "test_list.h"
#include "test.h"
#include "parameter_study.h"


NO_OPTIMIZATION

namespace testing{

using std::array;

void test_parameter_study()
{
    Tester tester = Tester("Test brute force parameter study utility");

    ADD_TEST(tester, "Test Range class",
        Const constant{1.0};
        ASSERT_EQUAL(constant[0], 1.0);
        ASSERT_EQUAL(constant[1], 1.0);
        ASSERT_EQUAL(constant[2], 1.0);

        LinearRange linrange{1, 16, 16};
        ASSERT_EQUAL(linrange[0], 1);
        ASSERT_EQUAL(linrange[15], 16);
        ASSERT_EQUAL(linrange[16], 16);
        ASSERT_EQUAL(linrange[17], 16);
        linrange = LinearRange{0, 1, 1};
        ASSERT_EQUAL(linrange[0], 0);
        ASSERT_EQUAL(linrange[1], 0);
        ASSERT_EQUAL(linrange[2], 0);
        linrange = LinearRange{0, 1, 0};
        ASSERT_EQUAL(linrange[0], 0);
        ASSERT_EQUAL(linrange[1], 0);
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 0);
        linrange = LinearRange{1, 0, 1};
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 1);
        ErrorHandler::clear_errors();

        PowRange powrange{1, 16, 16, 2};
        ASSERT_EQUAL(powrange[0], 1);
        ASSERT_EQUAL(powrange[15], 16);
        ASSERT_EQUAL(powrange[16], 16);
        ASSERT_EQUAL(powrange[17], 16);
        powrange = PowRange{0, 1, 1, 2};
        ASSERT_EQUAL(powrange[0], 0);
        ASSERT_EQUAL(powrange[1], 0);
        ASSERT_EQUAL(powrange[2], 0);
        powrange = PowRange{0, 1, 0, 2};
        ASSERT_EQUAL(powrange[0], 0);
        ASSERT_EQUAL(powrange[1], 0);

        ASSERT_TRUE((PowRange{1, 10, 5, 3.0}[1] < PowRange{1, 10, 5, 2.0}[1]));
        ASSERT_TRUE((PowRange{1, 10, 5, 3.0}[2] < PowRange{1, 10, 5, 2.0}[2]));
        ASSERT_TRUE((PowRange{1, 10, 5, 1.0}[1] < PowRange{1, 10, 5, 0.5}[1]));
    );

    ADD_TEST(tester, "Test ParameterStudy::get_total_combination_count()",
        ParameterStudy::Builder builder{
            .mechanism = Parameters::mechanism::chemkin_ar_he,
            .R_E = PowRange(1e-6, 10e-6, 5, 2),
            .rho_L = LinearRange(800.0, 1200.0, 5),
        };
        ParameterStudy ps1{builder};
        ASSERT_EQUAL(ps1.get_total_combination_count(), 5*5);

        builder.R_E = LinearRange(1e-6, 10e-6, 125);
        builder.rho_L = LinearRange(800.0, 1200.0, 125);
        builder.P_amb = LinearRange(101325.0, 101325.0, 1);
        builder.T_inf = LinearRange(293.15, 293.15, 125);
        builder.excitation_params[0] = LinearRange(-2.0e5, -2.0e5, 125);
        ParameterStudy ps2{builder};
        ASSERT_EQUAL(ps2.get_total_combination_count(), 125*125*125*125);
    );

    ADD_TEST(tester, "Test ParameterStudy::get_next_combination()",
        ParameterStudy::Builder builder{
            .mechanism = Parameters::mechanism::chemkin_ar_he,
            .R_E = LinearRange(0, 4, 5),
            .species = {"H2O"},
            .rho_L = LinearRange(0, 2, 3),
            .excitation_params = {
                LinearRange(0, 1, 2),
                Const(1),
                Const(1)
        }};
        ParameterStudy ps{builder};
        ASSERT_EQUAL(ps.get_total_combination_count(), 5*3*2);
        size_t ID = 0;
        for (size_t excitation_param = 0; excitation_param < 2; excitation_param++)
            for (size_t rho_L = 0; rho_L < 3; rho_L++)
                for (size_t R_E = 0; R_E < 5; R_E++)
                {
                    ASSERT_EQUAL(ps.get_next_combination_ID(), ID);
                    auto [success, cpar] = ps.get_next_combination();
                    ASSERT_TRUE(success);
                    ASSERT_EQUAL(cpar.ID, ID);
                    ASSERT_EQUAL(cpar.R_E, R_E);
                    ASSERT_EQUAL(cpar.rho_L, rho_L);
                    ASSERT_EQUAL(cpar.excitation_params[0], excitation_param);
                    ASSERT_EQUAL(cpar.excitation_params[1], 1);
                    ASSERT_EQUAL(cpar.excitation_params[2], 1);
                    ID++;
                }
        auto [success, cpar] = ps.get_next_combination();
        ASSERT_FALSE(success);
        ASSERT_EQUAL(cpar.error_ID, 0);
        ErrorHandler::clear_errors();
    );

    tester.run_tests();
}

}   // namespace testing

#endif // TEST
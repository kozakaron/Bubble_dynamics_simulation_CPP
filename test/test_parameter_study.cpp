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
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 0);
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

    ADD_TEST(tester, "Test ParameterCombinator::get_total_combination_count()",
        ParameterCombinator::Builder builder{
            .mechanism = Parameters::mechanism::chemkin_ar_he,
            .R_E = PowRange(1e-6, 10e-6, 5, 2),
            .rho_L = LinearRange(800.0, 1200.0, 5),
        };
        ParameterCombinator pc1{builder};
        ASSERT_EQUAL(pc1.get_total_combination_count(), 5*5);

        builder.R_E = LinearRange(1e-6, 10e-6, 125);
        builder.rho_L = LinearRange(800.0, 1200.0, 125);
        builder.P_amb = LinearRange(101325.0, 101325.0, 1);
        builder.T_inf = LinearRange(293.15, 293.15, 125);
        builder.excitation_params[0] = LinearRange(-2.0e5, -2.0e5, 125);
        ParameterCombinator pc2{builder};
        ASSERT_EQUAL(pc2.get_total_combination_count(), 125*125*125*125);
    );

    ADD_TEST(tester, "Test ParameterCombinator::get_next_combination()",
        ParameterCombinator::Builder builder{
            .mechanism = Parameters::mechanism::chemkin_ar_he,
            .R_E = LinearRange(0, 4, 5),
            .species = {"H2O"},
            .rho_L = LinearRange(0, 2, 3),
            .excitation_params = {
                LinearRange(0, 1, 2),
                Const(1),
                Const(1)
        }};
        ParameterCombinator pc{builder};
        ASSERT_EQUAL(pc.get_total_combination_count(), 5*3*2);
        size_t ID = 0;
        for (size_t excitation_param = 0; excitation_param < 2; excitation_param++)
            for (size_t rho_L = 0; rho_L < 3; rho_L++)
                for (size_t R_E = 0; R_E < 5; R_E++)
                {
                    ASSERT_EQUAL(pc.get_next_combination_ID(), ID);
                    auto [success, cpar] = pc.get_next_combination();
                    ASSERT_TRUE(success);
                    ASSERT_EQUAL(cpar.ID, ID);
                    ASSERT_EQUAL(cpar.R_E, R_E);
                    ASSERT_EQUAL(cpar.rho_L, rho_L);
                    ASSERT_EQUAL(cpar.excitation_params[0], excitation_param);
                    ASSERT_EQUAL(cpar.excitation_params[1], 1);
                    ASSERT_EQUAL(cpar.excitation_params[2], 1);
                    ID++;
                }
        auto [success, cpar] = pc.get_next_combination();
        ASSERT_FALSE(success);
        ASSERT_EQUAL(cpar.error_ID, ErrorHandler::no_error);
        ErrorHandler::clear_errors();
    );

    ADD_TEST(tester, "test SimulationData class",
        auto cpar = ControlParameters(ControlParameters::Builder{
            .ID                          = 0,
            .mechanism                   = Parameters::mechanism::chemkin_ar_he,
            .R_E                         = 1.00000000000000008e-05,    // bubble equilibrium radius [m]
            .species                     = {"O2"},
            .fractions                   = {1.00000000000000000e+00},
            .P_amb                       = 1.01325000000000000e+05,    // ambient pressure [Pa]
            .T_inf                       = 2.93149999999999977e+02,    // ambient temperature [K]
            .alfa_M                      = 3.49999999999999978e-01,    // water accommodation coefficient [-]
            .P_v                         = 2.33809999999999991e+03,    // vapour pressure [Pa]
            .mu_L                        = 1.00000000000000002e-03,    // dynamic viscosity [Pa*s]
            .rho_L                       = 9.98200000000000045e+02,    // liquid density [kg/m^3]
            .c_L                         = 1.48300000000000000e+03,    // sound speed [m/s]
            .surfactant                  = 1.00000000000000000e+00,    // surface tension modifier [-]
            .enable_heat_transfer        = true,
            .enable_evaporation          = true,
            .enable_reactions            = true,
            .enable_dissipated_energy    = true,
            .target_specie               = "H2O2",
            .excitation_params           = {-2.00000000000000000e+05, 3.00000000000000000e+04, 1.00000000000000000e+00},
            .excitation_type             = Parameters::excitation::sin_impulse
        });
    
        OdeSolution sol;
        sol.t = {0.0, 1.0};
        sol.x = {{
                    1.0000000000000001e-05, 0.0000000000000000e+00, 2.9314999999999998e+02, 0.0000000000000000e+00, 0.0000000000000000e+00,
                    0.0000000000000000e+00, 4.6517455752721046e-05, 0.0000000000000000e+00, 9.5926618412304955e-07, 0.0000000000000000e+00,
                    0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                    0.0000000000000000e+00
                }, {
                    1.0000850644975870e-05,  2.0736455319422995e-15,  2.9315000000000657e+02, -1.8915271945128421e-33,  2.4326921263456480e-10,
                   -1.9953268575847301e-22,  4.6494706600082013e-05,  1.0145497852167193e-25,  9.5926618412303917e-07,  0.0000000000000000e+00,
                    2.3380625551111345e-13,  2.2003345706478023e-08,  0.0000000000000000e+00,  0.0000000000000000e+00,  2.0732275271545515e-55,
                    1.1191854750923134e-06
                }};
        sol.num_dim = sol.x[0].size();

        SimulationData data(cpar, sol);
        ASSERT_APPROX(data.dissipated_energy, 1.1191854750923134e-06, 1e-20);
        ASSERT_APPROX(data.n_target_specie, 9.219091868668137e-17, 1e-5);
        ASSERT_APPROX(data.energy_demand, 356900.1798166697, 1e-5);

        ASSERT_EQUAL(SimulationData(cpar, OdeSolution()).energy_demand, SimulationData::infinite_energy_demand);
        const Parameters* par = Parameters::get_parameters(cpar.mechanism);
        cpar.target_specie = par->invalid_index;
        ASSERT_EQUAL(SimulationData(cpar, sol).energy_demand, SimulationData::infinite_energy_demand);
    );

    tester.run_tests();
}

}   // namespace testing

#endif // TEST
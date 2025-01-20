#ifdef TEST
#ifndef TEST_ODE_FUN_H
#define TEST_ODE_FUN_H
#include <cfloat>

#include "common.h"
#include "parameters.h"
#include "control_parameters.h"
#include "ode_fun.h"
#include "test.h"
#include "test_list.h"

NO_OPTIMIZATION

namespace testing{

using std::array;

class OdeFunTester_ar_he : public testing::Tester
{
public:
    ODE *ode;
    cpar_t* cpar;
    const Parameters *par;

    OdeFunTester_ar_he(std::string test_group_name): testing::Tester(test_group_name) {}

    void set_up() override
    {
        ErrorHandler::print_when_log = false;
        ErrorHandler::clear_errors();
        // Set up the ODE object
        ode = new ODE();
        cpar = new cpar_t();
        par = Parameters::get_parameters(Parameters::mechanism::chemkin_ar_he);

        cpar->ID = 0;
        cpar->mechanism = Parameters::mechanism::chemkin_ar_he;
        // Initial conditions:
        cpar->R_E = 10e-6;
        cpar->ratio = 1.0;
        cpar->set_species({par->get_species("O2")}, {1.0});
        // Ambient parameters:
        cpar->P_amb = 101325.0;
        cpar->T_inf = 293.15;
        // Liquid parameters:
        cpar->alfa_M = 0.35;
        cpar->P_v = 2338.1;
        cpar->mu_L = 0.001;
        cpar->rho_L = 998.2;
        cpar->c_L = 1483.0;
        cpar->surfactant = 1.0;
        // Simulation settings:
        cpar->enable_heat_transfer = true;
        cpar->enable_evaporation = true;
        cpar->enable_reactions = true;
        cpar->enable_dissipated_energy = true;
        cpar->target_specie = par->get_species("H2");
        // Excitation parameters:
        cpar->excitation_type = Parameters::excitation::sin_impulse;
        cpar->set_excitation_params({-2.0e5, 30000.0, 1.0});

        // Init the ODE object
        (void)ode->init(*cpar);
    }

    void tear_down() override
    {
        if (ErrorHandler::get_error_count() != 0) ErrorHandler::print_errors();
        ErrorHandler::clear_errors();
        delete ode;
        delete cpar;
    }
};

void test_ode_fun_ar_he()
{
    OdeFunTester_ar_he tester = OdeFunTester_ar_he("Test ode_fun_cpp.h's ODE class with chemkin_ar_he");

    const double t = 2.5255799280947053e-05;
    array<double, 12+4> x = {
        4.1866326662121154e-07, -1.0806071026205871e-02,  3.7345808358390459e+03,  2.6193400419650303e-05,  2.2272784074252325e-04,
        1.2216575395711056e-03,  6.1721264294022393e-01,  1.9348243023154722e-02,  2.9559164582684627e-01,  0.0000000000000000e+00,
        1.1445537293517848e-02,  5.5500381673767307e-03,  0.0000000000000000e+00,  0.0000000000000000e+00,  7.2195229829469272e-08,
        1.1381442517029756e-06
    };
    std::array<double, 12+4> dxdt;

    ADD_TEST(tester, "Test production_rate()",
        double T = 3734.580835839046;
        double M = 0.9506187582270826;
        double p = 29517684903.458748;
        array<double, 12> S = {
            1.6726028164481888e+09, 2.1116762033751166e+09, 2.1431788311718709e+09, 2.9341723434611363e+09, 2.6515733784052095e+09,
            2.9966686515624356e+09, 2.7504005615418749e+09, 3.5116915405108457e+09, 3.9333851262132778e+09, 2.0738884018602848e+09,
            1.7869604666669390e+09, 2.6488970129419308e+09
        };
        array<double, 12> H = {
            2.8942728456447578e+12,  1.1647431116371648e+12,  3.2122070817383286e+12,  1.2800798992004607e+12,  1.5473666390224338e+12,
            -7.1647582013218262e+11,  1.2008875056705388e+12,  1.9701027634066887e+12,  1.0646982960470887e+12,  7.1430166818375781e+11,
            7.1430166818375781e+11,  5.4183765876839316e+12
         };
        array<double, 7> M_eff_expected = {
            4.202460954083506 , 4.202460954083506 , 4.202460954083506 , 4.202460954083506 , 4.016609122669752 , 2.0327918805035803, 2.2060452245106026
        };
        array<double, 12> omega_dot_expected = {
            -1.3641074631640525e+05, -9.4837198306973488e+05, -3.6775705205366341e+06,  3.7956496431835137e+07, -4.6558168476294003e+07,
            4.9235108340222642e+07,  0.0000000000000000e+00, -2.5033468715366870e+07, -1.2422430271871166e+07,  0.0000000000000000e+00,
            0.0000000000000000e+00, -5.6423258621030436e+02
        };
        std::copy(S.begin(), S.end(), tester.ode->S);
        std::copy(H.begin(), H.end(), tester.ode->H);


        tester.ode->production_rate(T, M, p, x.data()+3);
        ASSERT_APPROX_ARRAY(tester.ode->M_eff, M_eff_expected, tester.par->num_third_bodies, 1e-15);
        ASSERT_APPROX_ARRAY(tester.ode->omega_dot, omega_dot_expected, tester.par->num_species, 1e-12);
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 0);
   
    );

    ADD_TEST(tester, "Test operator()",
        array<double, 12+4> dxdt_expected = {
            -1.0806071026205871e-02,  6.1535373020546820e+13, -1.6228795276785706e+12, -1.3640871809156740e+05, -9.4835473665913381e+05,
            -3.6774759243358579e+06,  3.8004288852543868e+07, -4.6556670290245622e+07, -3.3902018989090192e+08,  0.0000000000000000e+00,
            -2.5032582456895296e+07, -1.2422000517616974e+07,  0.0000000000000000e+00,  0.0000000000000000e+00, -5.6422699594115375e+02,
            6.9861519420936070e-04
        };
        tester.ode->cpar->enable_heat_transfer = true;
        tester.ode->cpar->enable_evaporation = true;
        tester.ode->cpar->enable_reactions = true;
        tester.ode->cpar->enable_dissipated_energy = true;

        const bool success = tester.ode->operator()(t, (const double*)x.data(), (double*)dxdt.data());
        ASSERT_APPROX_ARRAY(dxdt, dxdt_expected, tester.par->num_species+4, 1e-10);
        ASSERT_TRUE(success);
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 0);
    );

    tester.run_tests();
}

}   // namespace testing

#endif // TEST_ODE_FUN_H
#endif // TEST
#ifdef TEST
#ifndef TEST_PAR_CPAR_H
#define TEST_PAR_CPAR_H

#include "common.h"
#include "parameters.h"
#include "control_parameters.h"
#include "test.h"

NO_OPTIMIZATION

namespace testing{

void test_common()
{
    testing::Tester tester = testing::Tester("General tests");

    ADD_TEST(tester, "Test ErrorHandler class",
        ErrorHandler::print_when_log = false;
        const size_t idx1 = LOG_ERROR("Test error message", 33);
        const size_t idx2 = LOG_ERROR("Another test error message", 44);
        ASSERT_EQUAL(idx1, 0);
        ASSERT_EQUAL(idx2, 1);
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 2);
        ASSERT_EQUAL(ErrorHandler::get_error(0).ID, 33);
        ASSERT_EQUAL(ErrorHandler::get_error(1).ID, 44);
        ASSERT_TRUE(ErrorHandler::get_error(0).message == "Test error message");
        ErrorHandler::clear_errors();
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 0);
        ErrorHandler::print_when_log = true;
    );

    ADD_TEST(tester, "Test long double exists",
        long double x = 1.0L;
        ASSERT_TRUE(sizeof(x) > 8);    // long double is at least 80 bits
        long double fact_200 = 1.0L;
        for (int i = 1; i <= 200; ++i)
        {
            fact_200 *= i;
        }
        ASSERT_APPROX(fact_200, 7.88657867364790503383170119245e+374L, 1e-10);
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 0);
    );

    tester.run_tests();
}

void test_par_cpar()
{
    testing::Tester par_tester = testing::Tester("Test Parameters class");
    testing::Tester cpar_tester = testing::Tester("Test ControlParameters class");
    
    ADD_TEST(par_tester, "Test Parameters.get_parameters()",
        const Parameters* par;
        par = Parameters::get_parameters(Parameters::mechanism::chemkin_ar_he);
        ASSERT_TRUE(par->model == "chemkin_ar_he");
        ASSERT_EQUAL(par->num_elements, 6);
        ASSERT_EQUAL(par->num_species, 12);
        par = Parameters::get_parameters(Parameters::mechanism::chemkin_kaust2023_n2);
        ASSERT_TRUE(par->model == "chemkin_kaust2023_n2");
        ASSERT_EQUAL(par->num_elements, 6);
        ASSERT_EQUAL(par->num_species, 32);
        par = Parameters::get_parameters(Parameters::mechanism::chemkin_otomo2018_without_o);
        ASSERT_TRUE(par->model == "chemkin_otomo2018_without_o");
        ASSERT_EQUAL(par->num_elements, 5);
        ASSERT_EQUAL(par->num_species, 12);
        par = Parameters::get_parameters(Parameters::mechanism::chemkin_otomo2018);
        ASSERT_TRUE(par->model == "chemkin_otomo2018");
        ASSERT_EQUAL(par->num_elements, 5);
        ASSERT_EQUAL(par->num_species, 32);
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 0);

        ErrorHandler::print_when_log = false;
        par = Parameters::get_parameters((Parameters::mechanism)694269);
        ErrorHandler::print_when_log = true;
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 1);
        ErrorHandler::clear_errors();
    );

    ADD_TEST(par_tester, "Test getters",
        const Parameters* par;
        par = Parameters::get_parameters(Parameters::mechanism::chemkin_ar_he);
        ASSERT_TRUE("H2" == par->species_names[par->get_species("H2")]);
        ASSERT_TRUE("H2O" == par->species_names[par->get_species("H2O")]);
        ASSERT_TRUE("OH" == par->species_names[par->get_species("OH")]);
        ASSERT_TRUE("O2" == par->species_names[par->get_species("O2")]);
        (void)par->get_element("H");
        (void)par->get_element("O");
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 0);
    
        ErrorHandler::print_when_log = false;
        ASSERT_EQUAL(par->invalid_index, par->get_species("XX"));
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 1);
        ASSERT_EQUAL(par->invalid_index, par->get_element("X"))
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 2);
        ErrorHandler::print_when_log = true;
        ErrorHandler::clear_errors();
    );

    ADD_TEST(par_tester, "Access array elements",
        const Parameters* par;
        par = Parameters::get_parameters(Parameters::mechanism::chemkin_ar_he);
        (void)par->W[par->num_species-1];
        (void)par->lambdas[par->num_species-1];
        (void)par->temp_range[3*par->num_species-1];
        (void)par->a_low[par->num_species*(par->NASA_order+2)-1];
        (void)par->a_high[par->num_species*(par->NASA_order+2)-1];
        (void)par->A[par->num_reactions-1];
        (void)par->b[par->num_reactions-1];
        (void)par->E[par->num_reactions-1];
        ASSERT_APPROX(par->W[par->index_of_water], 18.01534, 1e-30);
    );


    ControlParameters cpar;
    const Parameters* par = Parameters::get_parameters(Parameters::mechanism::chemkin_otomo2018);

    cpar.ID = 0;
    cpar.mechanism = Parameters::mechanism::chemkin_otomo2018;
    // Initial conditions:
    cpar.R_E = 10e-6;
    cpar.ratio = 1.0;
    cpar.set_species({par->get_species("H2"), par->get_species("N2")}, {0.75, 0.25});
    // Ambient parameters:
    cpar.P_amb = 101325.0;
    cpar.T_inf = 293.15;
    // Liquid parameters:
    cpar.alfa_M = 0.35;
    cpar.P_v = 2338.1;
    cpar.mu_L = 0.001;
    cpar.rho_L = 998.2;
    cpar.c_L = 1483.0;
    cpar.surfactant = 1.0;
    // Simulation settings:
    cpar.enable_heat_transfer = true;
    cpar.enable_evaporation = true;
    cpar.enable_reactions = true;
    cpar.enable_dissipated_energy = true;
    cpar.target_specie = par->get_species("NH3");
    // Excitation parameters:
    cpar.excitation_type = Parameters::excitation::sin_impulse;
    cpar.set_excitation_params({-2.0e5, 30000.0, 1.0});

    ADD_TEST(cpar_tester, "Check preset cpar",
        ASSERT_TRUE(par->model == "chemkin_otomo2018");
        ASSERT_EQUAL(par->num_species, 32);
        ASSERT_TRUE(cpar.enable_heat_transfer);
        ASSERT_TRUE(cpar.enable_evaporation);
        ASSERT_TRUE(cpar.enable_reactions);
        ASSERT_TRUE(cpar.enable_dissipated_energy);
        ASSERT_EQUAL(cpar.target_specie, par->get_species("NH3"));
        ASSERT_EQUAL(cpar.excitation_type, Parameters::excitation::sin_impulse);
        ASSERT_EQUAL(cpar.excitation_params[0], -2.0e5);
        ASSERT_EQUAL(cpar.excitation_params[1], 30000.0);
        ASSERT_EQUAL(cpar.excitation_params[2], 1.0);
        ASSERT_EQUAL(cpar.species[0], par->get_species("H2"));
        ASSERT_EQUAL(cpar.species[1], par->get_species("N2"));
        ASSERT_EQUAL(cpar.fractions[0], 0.75);
        ASSERT_EQUAL(cpar.fractions[1], 0.25);
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 0);
    );

    ADD_TEST(cpar_tester, "Test cpar setters",
        // set_species()
        cpar.ID = 1;
        cpar.set_species({par->get_species("H2"), par->get_species("N2"), par->get_species("NH3")}, {0.5, 0.3, 0.2});
        ASSERT_EQUAL(cpar.n_species, 3);
        ASSERT_EQUAL(cpar.species[0], par->get_species("H2"));
        ASSERT_EQUAL(cpar.species[1], par->get_species("N2"));
        ASSERT_EQUAL(cpar.species[2], par->get_species("NH3"));
        ASSERT_EQUAL(cpar.fractions[0], 0.5);
        ASSERT_EQUAL(cpar.fractions[1], 0.3);
        ASSERT_EQUAL(cpar.fractions[2], 0.2);

        // set_excitation_params()
        cpar.excitation_type = Parameters::excitation::sin_impulse;
        cpar.set_excitation_params({-1.0e5, 20000.0, 0.5});
        ASSERT_EQUAL(cpar.excitation_params[0], -1.0e5);
        ASSERT_EQUAL(cpar.excitation_params[1], 20000.0);
        ASSERT_EQUAL(cpar.excitation_params[2], 0.5);

        // copy()
        cpar_t cpar_copy;
        cpar_copy.copy(cpar);
        ASSERT_EQUAL(cpar_copy.ID, 1);
        ASSERT_EQUAL(cpar_copy.n_species, 3);
        ASSERT_EQUAL(cpar_copy.species[0], par->get_species("H2"));
        ASSERT_EQUAL(cpar_copy.species[1], par->get_species("N2"));
        ASSERT_EQUAL(cpar_copy.species[2], par->get_species("NH3"));
        ASSERT_EQUAL(cpar_copy.fractions[0], 0.5);
        ASSERT_EQUAL(cpar_copy.fractions[1], 0.3);
        ASSERT_EQUAL(cpar_copy.fractions[2], 0.2);
        ASSERT_EQUAL(cpar_copy.excitation_type, Parameters::excitation::sin_impulse);
        ASSERT_EQUAL(cpar_copy.excitation_params[0], -1.0e5);
        ASSERT_EQUAL(cpar_copy.excitation_params[1], 20000.0);
        ASSERT_EQUAL(cpar_copy.excitation_params[2], 0.5);

        // error
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 0);
        cpar.set_species({par->get_species("H2"), par->get_species("N2")}, {0.75, 0.25, 0.0});
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 1);
        cpar.excitation_type = Parameters::excitation::no_excitation;
        cpar.set_excitation_params({-1.0e5, 20000.0, 0.5});
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 2);
        ErrorHandler::clear_errors();
    );

    ADD_TEST(cpar_tester, "Test cpar errors",
        ErrorHandler::clear_errors();
        cpar_t cpar;
        ASSERT_EQUAL(cpar.error_ID, ErrorHandler::no_error);
        cpar.set_species({par->get_species("H2"), par->get_species("N2")}, {0.75, 0.255});
        ASSERT_FALSE(cpar.error_ID == ErrorHandler::no_error);
        ASSERT_EQUAL(cpar.error_ID, 0);
        Error error = ErrorHandler::get_error(cpar.error_ID);
        ASSERT_TRUE(error.message.find("1.0") != std::string::npos);
        ErrorHandler::clear_errors();
    );


    par_tester.run_tests();
    cpar_tester.run_tests();
}

}   // namespace testing

#endif  // TEST_PAR_CPAR_H
#endif  // TEST
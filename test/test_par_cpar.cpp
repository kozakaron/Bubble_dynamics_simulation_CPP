#ifdef TEST

#include "nlohmann/json.hpp"
#include "common.h"
#include "parameters.h"
#include "control_parameters.h"
#include "test.h"


using ordered_json = nlohmann::ordered_json;

NO_OPTIMIZATION

namespace testing{

void test_common()
{
    testing::Tester tester = testing::Tester("General tests");

    ADD_TEST(tester, "Test ErrorHandler class",
        ErrorHandler::print_when_log = false;
        const size_t idx1 = LOG_ERROR("Test error message");
        const size_t idx2 = LOG_ERROR("Another test error message", 44);
        const size_t idx3 = LOG_ERROR(Error::severity::warning, Error::type::preprocess, "Another test error message");
        const size_t idx4 = LOG_ERROR(Error::severity::info, Error::type::postprocess, "Another test error message", 66);
        ASSERT_EQUAL(idx1, 0);
        ASSERT_EQUAL(idx2, 1);
        ASSERT_EQUAL(idx3, 2);
        ASSERT_EQUAL(idx4, 3);
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 4);
        ASSERT_EQUAL(ErrorHandler::get_error(0).ID, 0);
        ASSERT_EQUAL(ErrorHandler::get_error(1).ID, 44);
        ASSERT_EQUAL(ErrorHandler::get_error(2).ID, 0);
        ASSERT_EQUAL(ErrorHandler::get_error(3).ID, 66);
        ASSERT_TRUE(ErrorHandler::get_error(0).message == "Test error message");
        ASSERT_TRUE(ErrorHandler::get_error(1).message == "Another test error message");
        ASSERT_EQUAL(ErrorHandler::get_error(0).error_severity, Error::severity::error);
        ASSERT_EQUAL(ErrorHandler::get_error(1).error_severity, Error::severity::error);
        ASSERT_EQUAL(ErrorHandler::get_error(2).error_severity, Error::severity::warning);
        ASSERT_EQUAL(ErrorHandler::get_error(3).error_severity, Error::severity::info);
        ASSERT_EQUAL(ErrorHandler::get_error(0).error_type, Error::type::general);
        ASSERT_EQUAL(ErrorHandler::get_error(1).error_type, Error::type::general);
        ASSERT_EQUAL(ErrorHandler::get_error(2).error_type, Error::type::preprocess);
        ASSERT_EQUAL(ErrorHandler::get_error(3).error_type, Error::type::postprocess);
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

    const Parameters* par = Parameters::get_parameters(Parameters::mechanism::chemkin_otomo2018);
    ControlParameters cpar{ControlParameters::Builder{
        .ID                          = 0,
        .mechanism                   = Parameters::mechanism::chemkin_otomo2018,
        .R_E                         = 1.00000000000000008e-05,    // bubble equilibrium radius [m]
        .species                     = {"H2", "N2"},
        .fractions                   = {7.50000000000000000e-01, 2.50000000000000000e-01},
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
        .target_specie               = "NH3",
        .excitation_params           = {-2.00000000000000000e+05, 3.00000000000000000e+04, 1.00000000000000000e+00},
        .excitation_type             = Parameters::excitation::sin_impulse
    }};

    ADD_TEST(cpar_tester, "Check preset cpar",
        ASSERT_TRUE(par->model == "chemkin_otomo2018");
        ASSERT_EQUAL(par->num_species, 32);
        ASSERT_EQUAL(cpar.mechanism, Parameters::mechanism::chemkin_otomo2018);
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
        ASSERT_EQUAL(cpar.num_initial_species, 3);
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

        // copy
        ControlParameters cpar_copy;
        cpar_copy = cpar;
        ASSERT_EQUAL(cpar_copy.ID, 1);
        ASSERT_EQUAL(cpar_copy.num_initial_species, 3);
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
        ControlParameters cpar{ControlParameters::Builder{
            .ID                          = 69,
            .mechanism                   = Parameters::mechanism::chemkin_otomo2018,
        }};
        ASSERT_EQUAL(cpar.error_ID, ErrorHandler::no_error);
        cpar.set_species({par->get_species("H2"), par->get_species("N2")}, {0.75, 0.255});
        ASSERT_FALSE(cpar.error_ID == ErrorHandler::no_error);
        ASSERT_EQUAL(cpar.error_ID, 0);
        Error error = ErrorHandler::get_error(cpar.error_ID);
        ASSERT_TRUE(error.message.find("1.0") != std::string::npos);
        ASSERT_EQUAL(error.ID, 69);
        ErrorHandler::clear_errors();
    );

    ADD_TEST(cpar_tester, "Test cpar JSON-based constructor",
        ErrorHandler::clear_errors();
    
        // Valid JSON input (no default values, arbitrary data)
        ordered_json valid_json = {
            {"ID", 123},
            {"mechanism", "Chemkin_Otomo2018"}, // Case insensitive
            {"R_E", 0.12345},
            {"species", {"H2", "O2"}},
            {"fractions", {0.6, 0.4}},
            {"P_amb", 50000.0},
            {"T_inf", 350.0},
            {"alfa_M", 0.25},
            {"P_v", 1500.0},
            {"mu_L", 0.002},
            {"rho_L", 800.0},
            {"c_L", 1200.0},
            {"surfactant", 0.8},
            {"enable_heat_transfer", false},
            {"enable_evaporation", true},
            {"enable_reactions", false},
            {"enable_dissipated_energy", true},
            {"target_specie", "h2o"},   // Case insensitive
            {"excitation_params", {-100000.0, 25000.0, 0.75}},
            {"excitation_type", "Sin_Impulse"} // Case insensitive
        };
    
        ControlParameters cpar(valid_json);
        ASSERT_EQUAL(cpar.ID, 123);
        ASSERT_EQUAL(cpar.mechanism, Parameters::mechanism::chemkin_otomo2018);
        ASSERT_EQUAL(cpar.R_E, 0.12345);
        ASSERT_EQUAL(cpar.species[0], par->get_species("H2"));
        ASSERT_EQUAL(cpar.species[1], par->get_species("O2"));
        ASSERT_EQUAL(cpar.fractions[0], 0.6);
        ASSERT_EQUAL(cpar.fractions[1], 0.4);
        ASSERT_EQUAL(cpar.P_amb, 50000.0);
        ASSERT_EQUAL(cpar.T_inf, 350.0);
        ASSERT_EQUAL(cpar.alfa_M, 0.25);
        ASSERT_EQUAL(cpar.P_v, 1500.0);
        ASSERT_EQUAL(cpar.mu_L, 0.002);
        ASSERT_EQUAL(cpar.rho_L, 800.0);
        ASSERT_EQUAL(cpar.c_L, 1200.0);
        ASSERT_EQUAL(cpar.surfactant, 0.8);
        ASSERT_FALSE(cpar.enable_heat_transfer);
        ASSERT_TRUE(cpar.enable_evaporation);
        ASSERT_FALSE(cpar.enable_reactions);
        ASSERT_TRUE(cpar.enable_dissipated_energy);
        ASSERT_EQUAL(cpar.target_specie, par->get_species("H2O"));
        ASSERT_EQUAL(cpar.excitation_type, Parameters::excitation::sin_impulse);
        ASSERT_EQUAL(cpar.excitation_params[0], -100000.0);
        ASSERT_EQUAL(cpar.excitation_params[1], 25000.0);
        ASSERT_EQUAL(cpar.excitation_params[2], 0.75);
    
        // Missing values
        ordered_json missing_values_json = {
            {"ID", 124},
            {"mechanism", "chemkin_ar_he"} // Missing other fields
        };
    
        ControlParameters cpar_missing(missing_values_json);
        ControlParameters default_cpar;
        ASSERT_EQUAL(cpar_missing.ID, 124);
        ASSERT_EQUAL(cpar_missing.mechanism, Parameters::mechanism::chemkin_ar_he);
        ASSERT_EQUAL(cpar_missing.R_E, default_cpar.R_E);
        ASSERT_EQUAL(cpar_missing.num_initial_species, default_cpar.num_initial_species);
        ASSERT_EQUAL(cpar_missing.P_amb, default_cpar.P_amb);
        ASSERT_EQUAL(cpar_missing.T_inf, default_cpar.T_inf);
        ASSERT_EQUAL(cpar_missing.alfa_M, default_cpar.alfa_M);
        ASSERT_EQUAL(cpar_missing.excitation_type, default_cpar.excitation_type);
    
        // Invalid mechanism
        ordered_json invalid_mechanism_json = {
            {"ID", 125},
            {"mechanism", "Invalid_Mechanism"}
        };
    
        ErrorHandler::print_when_log = false;
        ControlParameters cpar_invalid_mechanism(invalid_mechanism_json);
        ASSERT_EQUAL(cpar_invalid_mechanism.mechanism, default_cpar.mechanism);
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 1);
        ASSERT_TRUE(ErrorHandler::get_error(0).message.find("Invalid mechanism") != std::string::npos);
    
        // Invalid excitation type
        ordered_json invalid_excitation_json = {
            {"ID", 126},
            {"excitation_type", "Invalid_Excitation"}
        };
    
        ControlParameters cpar_invalid_excitation(invalid_excitation_json);
        ASSERT_EQUAL(cpar_invalid_excitation.excitation_type, default_cpar.excitation_type);
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 2);
        ASSERT_TRUE(ErrorHandler::get_error(1).message.find("Invalid excitation type") != std::string::npos);
    
        ErrorHandler::print_when_log = true;
        ErrorHandler::clear_errors();
    );


    par_tester.run_tests();
    cpar_tester.run_tests();
}

}   // namespace testing

#endif  // TEST
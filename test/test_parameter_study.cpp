#ifdef TEST
#include <array>

#include "common.h"
#include "test_list.h"
#include "test.h"
#include "parameters.h"
#include "parameter_combinator.h"
#include "parameter_study.h"
#include "nlohmann/json.hpp"

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
        ASSERT_EQUAL(linrange[0], 1.0);
        ASSERT_EQUAL(linrange[15], 16.0);
        ASSERT_EQUAL(linrange[16], 16.0);
        ASSERT_EQUAL(linrange[17], 16.0);
        linrange = LinearRange{0, 1, 1};
        ASSERT_EQUAL(linrange[0], 0.0);
        ASSERT_EQUAL(linrange[1], 0.0);
        ASSERT_EQUAL(linrange[2], 0.0);
        linrange = LinearRange{0, 1, 0};
        ASSERT_EQUAL(linrange[0], 0.0);
        ASSERT_EQUAL(linrange[1], 0.0);
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 0);
        linrange = LinearRange{1, 0, 1};
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 0);
        ErrorHandler::clear_errors();

        LogRange logrange{1, 16, 16};
        ASSERT_NEAR(logrange[0], 1.0, 1e-10);
        ASSERT_NEAR(logrange[15], 16.0, 1e-10);
        ASSERT_NEAR(logrange[16], 16.0, 1e-10);
        ASSERT_NEAR(logrange[17], 16.0, 1e-10);
        logrange = LogRange{1, 2, 1};
        ASSERT_NEAR(logrange[0], 1.0, 1e-10);
        ASSERT_NEAR(logrange[1], 1.0, 1e-10);
        ASSERT_NEAR(logrange[2], 1.0, 1e-10);
        logrange = LogRange{1, 2, 0};
        ASSERT_NEAR(logrange[0], 1.0, 1e-10);
        ASSERT_NEAR(logrange[1], 1.0, 1e-10);

        // values from Python's numpy.logspace(1, 10, 5)
        logrange = LogRange{10, 1e10, 5};
        ASSERT_APPROX(logrange[0], 1.00000000e+01, 1e-5);
        ASSERT_APPROX(logrange[1], 1.77827941e+03, 1e-5);
        ASSERT_APPROX(logrange[2], 3.16227766e+05, 1e-5);
        ASSERT_APPROX(logrange[3], 5.62341325e+07, 1e-5);
        ASSERT_APPROX(logrange[4], 1.00000000e+10, 1e-5);

        logrange = LogRange{-10, -1e10, 5};
        ASSERT_APPROX(logrange[0], -1.00000000e+01, 1e-5);
        ASSERT_APPROX(logrange[1], -1.77827941e+03, 1e-5);
        ASSERT_APPROX(logrange[2], -3.16227766e+05, 1e-5);
        ASSERT_APPROX(logrange[3], -5.62341325e+07, 1e-5);
        ASSERT_APPROX(logrange[4], -1.00000000e+10, 1e-5);

        logrange = LogRange{0, 1, 5};
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 1);
        logrange = LogRange{-1, 1, 5};
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 2);
        logrange = LogRange{1, -1, 5};
        ASSERT_EQUAL(ErrorHandler::get_error_count(), 3);
        ErrorHandler::clear_errors();

        GeomRange geomrange{-1, 16, 16, 2.0};
        ASSERT_NEAR(geomrange[0], -1.0, 1e-10);
        ASSERT_NEAR(geomrange[15], 16.0, 1e-10);
        ASSERT_NEAR(geomrange[16], 16.0, 1e-10);
        ASSERT_NEAR(geomrange[17], 16.0, 1e-10);
        double a1 = geomrange[1] - geomrange[0];
        double a2 = geomrange[2] - geomrange[1];
        double a3 = geomrange[3] - geomrange[2];
        ASSERT_TRUE(a3 > a2 && a2 > a1);
        ASSERT_NEAR(a2/a1, 2.0, 1e-10);
        ASSERT_NEAR(a3/a2, 2.0, 1e-10);

        geomrange = GeomRange{1, -16, 16, 0.5};
        ASSERT_NEAR(geomrange[0], 1.0, 1e-10);
        ASSERT_NEAR(geomrange[15], -16.0, 1e-10);
        ASSERT_NEAR(geomrange[16], -16.0, 1e-10);
        ASSERT_NEAR(geomrange[17], -16.0, 1e-10);
        a1 = std::abs(geomrange[1] - geomrange[0]);
        a2 = std::abs(geomrange[2] - geomrange[1]);
        a3 = std::abs(geomrange[3] - geomrange[2]);
        ASSERT_TRUE(a3 < a2 && a2 < a1);
        ASSERT_NEAR(a2/a1, 0.5, 1e-10);
        ASSERT_NEAR(a3/a2, 0.5, 1e-10);
    );

    ADD_TEST(tester, "Test Range class JSON serialization",
        // Test Const JSON conversion
        Const const_range{5.0};
        auto const_json = const_range.to_json();
        ASSERT_TRUE(const_json.at("type").get<std::string>() == "Const");
        ASSERT_EQUAL(const_json.at("value").get<double>(), 5.0);

        // Test LinearRange JSON conversion
        LinearRange linear_range{1.0, 10.0, 5};
        auto linear_json = linear_range.to_json();
        ASSERT_TRUE(linear_json.at("type").get<std::string>() == "LinearRange");
        ASSERT_EQUAL(linear_json.at("start").get<double>(), 1.0);
        ASSERT_EQUAL(linear_json.at("end").get<double>(), 10.0);
        ASSERT_EQUAL(linear_json.at("num_steps").get<int>(), 5);

        // Test LogRange JSON conversion
        LogRange log_range{1.0, 100.0, 3};
        auto log_json = log_range.to_json();
        ASSERT_TRUE(log_json.at("type").get<std::string>() == "LogRange");
        ASSERT_EQUAL(log_json.at("start").get<double>(), 1.0);
        ASSERT_EQUAL(log_json.at("end").get<double>(), 100.0);
        ASSERT_EQUAL(log_json.at("num_steps").get<int>(), 3);

        // Test GeomRange JSON conversion
        GeomRange geom_range{1.0, 16.0, 5, 2.0};
        auto geom_json = geom_range.to_json();
        ASSERT_TRUE(geom_json.at("type").get<std::string>() == "GeomRange");
        ASSERT_EQUAL(geom_json.at("start").get<double>(), 1.0);
        ASSERT_EQUAL(geom_json.at("end").get<double>(), 16.0);
        ASSERT_EQUAL(geom_json.at("num_steps").get<int>(), 5);
        ASSERT_EQUAL(geom_json.at("q").get<double>(), 2.0);
    );

    ADD_TEST(tester, "Test ParameterCombinator::get_total_combination_count()",
        ParameterCombinator::Builder builder{
            .mechanism = Parameters::mechanism::chemkin_ar_he,
            .R_E = LogRange(1e-6, 10e-6, 5),
            .rho_L = LinearRange(800.0, 1200.0, 5),
        };
        ParameterCombinator pc1{builder};
        ASSERT_EQUAL(pc1.get_total_combination_count(), 5*5);

        builder.R_E = LinearRange(1e-6, 10e-6, 125);
        builder.rho_L = LinearRange(800.0, 1200.0, 125);
        builder.P_amb = LinearRange(101325.0, 101325.0, 1);
        builder.T_inf = LinearRange(293.15, 293.15, 125);
        builder.excitation_params[0] = GeomRange(-2.0e5, -2.0e5, 125, 1.0);
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
        for (size_t rho_L = 0; rho_L < 3; rho_L++)
            for (size_t R_E = 0; R_E < 5; R_E++)
                for (size_t excitation_param = 0; excitation_param < 2; excitation_param++)
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

    ADD_TEST(tester, "Test ParameterCombinator JSON constructor",
        std::string json_str = R"({
            "mechanism": "chemkin_ar_he",
            "R_E": {"type": "Linear Range", "start": 1e-6, "end": 10e-6, "num_steps": 3},
            "species": ["O2", "H2"],
            "fractions": [0.7, 0.3],
            "P_amb": {"type": "Const", "value": 101325.0},
            "T_inf": {"type": "logRange", "start": 250.0, "end": 350.0, "num_steps": 4},
            "alfa_M": {"type": "GeomRange", "start": 0.1, "end": 0.9, "num_steps": 5, "q": 1.5},
            "P_v": {"type": "CONST", "value": 2338.1},
            "mu_L": {"type": "Const", "value": 0.001},
            "rho_L": {"type": "LinearRange", "start": 900.0, "end": 1100.0, "num_steps": 2},
            "c_L": {"type": "Const", "value": 1483.0},
            "surfactant": {"type": "Const", "value": 1.0},
            "enable_heat_transfer": true,
            "enable_evaporation": false,
            "enable_reactions": true,
            "enable_dissipated_energy": false,
            "target_specie": "H2O2",
            "excitation_params": [
                {"type": "Const", "value": -200000.0},
                {"type": "LinearRange", "start": 20000.0, "end": 40000.0, "num_steps": 2},
                {"type": "Const", "value": 1.0}
            ],
            "excitation_type": "sin_impulse"
        })";

        auto json_obj = nlohmann::ordered_json::parse(json_str);
        ParameterCombinator pc_from_json(json_obj);
        ASSERT_EQUAL(pc_from_json.get_total_combination_count(), 3 * 4 * 5 * 2 * 2);
        
        auto [success1, cpar1] = pc_from_json.get_next_combination();
        ASSERT_TRUE(success1);
        ASSERT_EQUAL(cpar1.ID, 0);
        ASSERT_EQUAL(cpar1.mechanism, Parameters::mechanism::chemkin_ar_he);
        ASSERT_EQUAL(cpar1.num_initial_species, 2);
        const Parameters* par = Parameters::get_parameters(Parameters::mechanism::chemkin_ar_he);
        ASSERT_EQUAL(cpar1.species[0], par->get_species("O2"));
        ASSERT_EQUAL(cpar1.species[1], par->get_species("H2"));
        ASSERT_EQUAL(cpar1.fractions[0], 0.7);
        ASSERT_EQUAL(cpar1.fractions[1], 0.3);
        ASSERT_EQUAL(cpar1.P_amb, 101325.0);
        ASSERT_EQUAL(cpar1.enable_heat_transfer, true);
        ASSERT_EQUAL(cpar1.enable_evaporation, false);
        ASSERT_EQUAL(cpar1.enable_reactions, true);
        ASSERT_EQUAL(cpar1.enable_dissipated_energy, false);
        ASSERT_EQUAL(cpar1.target_specie, par->get_species("H2O2"));
        ASSERT_EQUAL(cpar1.excitation_type, Parameters::excitation::sin_impulse);
        ASSERT_EQUAL(cpar1.excitation_params[0], -200000.0);
        ASSERT_EQUAL(cpar1.excitation_params[1], 20000.0);  // First value from LinearRange
        ASSERT_EQUAL(cpar1.excitation_params[2], 1.0);
        
        auto [success2, cpar2] = pc_from_json.get_next_combination();
        ASSERT_TRUE(success2);
        ASSERT_EQUAL(cpar2.ID, 1);
        ASSERT_EQUAL(cpar2.excitation_params[1], 40000.0);  // Second value from LinearRange
        
        // Test JSON serialization
        auto serialized_json = pc_from_json.to_json();
        ASSERT_TRUE(serialized_json.contains("parameter_study"));
        ASSERT_TRUE(serialized_json["parameter_study"]["mechanism"].get<std::string>() == "chemkin_ar_he");
        ASSERT_EQUAL(serialized_json["parameter_study"]["species"].size(), 2);
        ASSERT_TRUE(serialized_json["parameter_study"]["target_specie"].get<std::string>() == "H2O2");

        // Test error handling - invalid JSON structure
        ErrorHandler::clear_errors();
        std::string invalid_json = R"({"invalid": "structure"})";
        auto invalid_json_obj = nlohmann::ordered_json::parse(invalid_json);
        ParameterCombinator pc_invalid(invalid_json_obj);
        ASSERT_TRUE(ErrorHandler::get_error_count() > 0);
        ErrorHandler::clear_errors();

        // Test error handling - test invalid range type
        auto invalid_range_obj = json_obj;
        invalid_range_obj["R_E"] = nlohmann::json{{{"type", "InvalidRange"}}, {"value", 1e-6}};
        ParameterCombinator pc_invalid_range(invalid_range_obj);
        ASSERT_TRUE(ErrorHandler::get_error_count() > 0);
        ErrorHandler::clear_errors();
        
        // Has const default values
        std::string minimal_json = R"({
            "mechanism": "chemkin_ar_he"
        })";
        ErrorHandler::print_when_log = false;
        auto minimal_json_obj = nlohmann::ordered_json::parse(minimal_json);
        ParameterCombinator pc_minimal(minimal_json_obj);
        ASSERT_EQUAL(pc_minimal.get_total_combination_count(), 1);  // All defaults are Const
        ErrorHandler::print_when_log = true;
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
    
        SimulationData data(cpar);
        data.sol.t = {0.0, 1.0};
        data.sol.x = {{
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
        data.sol.num_dim = data.sol.x[0].size();
        data.postprocess();

        ASSERT_APPROX(data.dissipated_energy, 1.1191854750923134e-06, 1e-10);
        ASSERT_APPROX(data.n_target_specie, 9.219091868668137e-17, 1e-5);
        ASSERT_APPROX(data.energy_demand, 356900.1798166697, 1e-5);

        ASSERT_EQUAL(SimulationData(cpar).energy_demand, SimulationData::infinite_energy_demand);
    );

    tester.run_tests();
}

}   // namespace testing

#endif // TEST
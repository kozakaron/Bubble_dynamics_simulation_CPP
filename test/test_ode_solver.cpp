#ifndef TEST
#define TEST // TODO: remove this line when the test is complete
#endif

#ifdef TEST
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

#include "common.h"
#include "test_list.h"
#include "test.h"
#include "ode_solver.h"
#include "ode_solver_sundials.h"
#include "parameter_study.h"


NO_OPTIMIZATION
namespace testing{
using std::vector;
using std::string;


struct TestReferenceData
{
    const string name;
    const ControlParameters cpar;
    const double t_final;
    const double T_final;
    const vector<double> mol_fractions_final;
};

vector<TestReferenceData> init_test_data();
void test_CVODE_solver(vector<TestReferenceData> &test_data);

void test_ode_solver()
{
    vector<TestReferenceData> test_data = init_test_data();
    test_CVODE_solver(test_data);
}


vector<TestReferenceData> init_test_data()
{
    vector<TestReferenceData> test_data;
    test_data.push_back(TestReferenceData{
        .name = "Test chemkin_kaust2023_n2; 1000 bar; sin_impulse; ",
        .cpar = ControlParameters{ControlParameters::Builder{
            .mechanism                               = Parameters::mechanism::chemkin_kaust2023_n2,
            .R_E                                     = 50.0e-6,
            .species                                 = {"H2", "O2", "N2", "AR", "HE"},
            .fractions                               = {0.3, 0.1, 0.15, 0.2, 0.25},
            .P_amb                                   = 1000 * 100000.0,
            .T_inf                                   = 2293.15,
            .alfa_M                                  = 0.35,
            .P_v                                     = 2338.339978450019,
            .mu_L                                    = 0.001,
            .rho_L                                   = 998.2,
            .c_L                                     = 1483.0,
            .surfactant                              = 1.0,
            .enable_heat_transfer                    = false,
            .enable_evaporation                      = false,
            .enable_reactions                        = true,
            .enable_dissipated_energy                = false,
            .target_specie                           = "NH3",
            .excitation_params                       = {-0.3e5, 25000.0, 1.0},
            .excitation_type                         = Parameters::excitation::sin_impulse
        }},
        .t_final = 0.00002,
        .T_final = 3708.918,
        .mol_fractions_final = {
            2.1968250000000000e-01,    // AR
            9.8731380000000001e-03,    // H
            1.1256170000000000e-01,    // H2
            9.0446219999999997e-10,    // H2NN
            2.6828989999999998e-07,    // H2NO
            2.0696580000000001e-01,    // H2O
            4.8871299999999997e-06,    // H2O2
            2.7460309999999999e-01,    // HE
            8.9944809999999997e-06,    // HNO
            9.9938880000000007e-09,    // HNO2
            8.8758649999999995e-08,    // HNOH
            7.6383519999999999e-06,    // HO2
            2.6810990000000001e-08,    // HON
            3.3442499999999998e-07,    // HONO
            6.2639469999999994e-11,    // HONO2
            7.1758750000000001e-06,    // N
            1.6381430000000000e-01,    // N2
            1.7192440000000001e-08,    // N2H2
            3.9088070000000001e-10,    // N2H3
            7.8870879999999999e-11,    // N2H4
            1.3344939999999999e-06,    // N2O
            1.3187289999999999e-05,    // NH
            2.4648820000000002e-05,    // NH2
            6.7023719999999995e-08,    // NH2OH
            7.2750680000000000e-05,    // NH3
            5.2599770000000000e-06,    // NNH
            1.7536150000000000e-03,    // NO
            7.1548869999999995e-07,    // NO2
            4.8678139999999999e-12,    // NO3
            4.6630880000000002e-04,    // O
            3.2662369999999997e-04,    // O2
            9.8055799999999995e-03     // OH
        }
    });

    test_data.push_back(TestReferenceData{
        .name = "Test chemkin_otomo2018; 1000 bar; sin_impulse; ",
        .cpar = ControlParameters{ControlParameters::Builder{
            .mechanism                               = Parameters::mechanism::chemkin_otomo2018,
            .R_E                                     = 50.0e-6,
            .species                                 = {"H2", "O2", "N2", "AR", "HE"},
            .fractions                               = {0.3, 0.1, 0.15, 0.2, 0.25},
            .P_amb                                   = 1000 * 100000.0,
            .T_inf                                   = 2293.15,
            .alfa_M                                  = 0.35,
            .P_v                                     = 2338.339978450019,
            .mu_L                                    = 0.001,
            .rho_L                                   = 998.2,
            .c_L                                     = 1483.0,
            .surfactant                              = 1.0,
            .enable_heat_transfer                    = false,
            .enable_evaporation                      = false,
            .enable_reactions                        = true,
            .enable_dissipated_energy                = false,
            .target_specie                           = "NH3",
            .excitation_params                       = {-0.3e5, 25000.0, 1.0},
            .excitation_type                         = Parameters::excitation::sin_impulse
        }},
        .t_final = 0.000019,
        .T_final = 3708.459,
        .mol_fractions_final = {
            1.7147530000000001e-03,    // NO
            4.7374340000000001e-05,    // NH3
            1.1261320000000000e-01,    // H2
            3.2553060000000002e-04,    // O2
            9.8612890000000005e-03,    // H
            4.6530829999999999e-04,    // O
            9.9167700000000001e-03,    // OH
            7.6818309999999997e-06,    // HO2
            2.0688400000000001e-01,    // H2O
            4.9154240000000002e-06,    // H2O2
            2.5173109999999999e-05,    // NH2
            1.3538499999999999e-05,    // NH
            7.2717080000000004e-06,    // N
            5.6963080000000002e-06,    // NNH
            6.7092530000000006e-08,    // NH2OH
            1.4312049999999999e-07,    // H2NO
            8.8761020000000003e-08,    // HNOH
            8.7074390000000007e-06,    // HNO
            4.4122959999999998e-08,    // HON
            8.3245539999999997e-07,    // NO2
            6.6577550000000000e-07,    // HONO
            9.9831850000000002e-09,    // HNO2
            4.8556990000000002e-12,    // NO3
            6.2515369999999997e-11,    // HONO2
            1.4659150000000000e-06,    // N2O
            1.7164420000000001e-10,    // N2H4
            1.8342140000000001e-09,    // N2H3
            2.0136900000000001e-08,    // N2H2
            2.0197099999999998e-09,    // H2NN
            2.1967059999999999e-01,    // AR
            2.7458830000000001e-01,    // HE
            1.6383650000000000e-01,    // N2
        }
    });

    test_data.push_back(TestReferenceData{
        .name = "Test chemkin_otomo2018; 1 bar; sin_impulse; ",
        .cpar = ControlParameters{ControlParameters::Builder{
            .mechanism                               = Parameters::mechanism::chemkin_otomo2018,
            .R_E                                     = 50.0e-6,
            .species                                 = {"H2", "O2", "N2", "AR", "HE"},
            .fractions                               = {0.3, 0.1, 0.15, 0.2, 0.25},
            .P_amb                                   = 1 * 100000.0,
            .T_inf                                   = 2293.15,
            .alfa_M                                  = 0.35,
            .P_v                                     = 2338.339978450019,
            .mu_L                                    = 0.001,
            .rho_L                                   = 998.2,
            .c_L                                     = 1483.0,
            .surfactant                              = 1.0,
            .enable_heat_transfer                    = false,
            .enable_evaporation                      = false,
            .enable_reactions                        = true,
            .enable_dissipated_energy                = false,
            .target_specie                           = "NH3",
            .excitation_params                       = {-0.3e5, 25000.0, 1.0},
            .excitation_type                         = Parameters::excitation::sin_impulse
        }},
        .t_final = 0.00045,
        .T_final = 3055.213,
        .mol_fractions_final = {
            2.1745839999999998e-03,    // NO
            4.9901519999999998e-08,    // NH3
            1.1182560000000000e-01,    // H2
            5.0010740000000003e-03,    // O2
            5.8488129999999999e-02,    // H
            9.0319649999999994e-03,    // O
            3.0766260000000000e-02,    // OH
            3.6937060000000001e-06,    // HO2
            1.5696050000000000e-01,    // H2O
            2.0869440000000000e-07,    // H2O2
            1.5075150000000001e-07,    // NH2
            7.1080589999999995e-07,    // NH
            4.8429830000000002e-06,    // N
            2.9867929999999997e-08,    // NNH
            7.3591389999999999e-12,    // NH2OH
            2.5580699999999999e-10,    // H2NO
            1.1977389999999999e-10,    // HNOH
            2.9965900000000000e-07,    // HNO
            1.5676690000000000e-12,    // HON
            1.9966610000000000e-07,    // NO2
            1.0753880000000000e-08,    // HONO
            1.3011310000000001e-10,    // HNO2
            9.9516089999999999e-14,    // NO3
            1.8123980000000000e-13,    // HONO2
            7.7678420000000004e-08,    // N2O
            4.6763430000000000e-17,    // N2H4
            1.6458940000000002e-14,    // N2H3
            4.1135240000000003e-12,    // N2H2
            2.5899520000000001e-13,    // H2NN
            2.0894409999999999e-01,    // AR
            2.6118010000000003e-01,    // HE
            1.5561749999999999e-01,    // N2
        }
    });


    for (auto &data : test_data)
    {
        double sum_mol_fraction = 0.0;
        for (auto &mol_fraction : data.mol_fractions_final)
        {
            sum_mol_fraction += mol_fraction;
        }
        if (std::abs(1.0 - sum_mol_fraction) > 1.0e-6)
        {
            LOG_ERROR("Sum of molar fractions is not 1.0 but " + std::to_string(sum_mol_fraction));
        }
    }

    return test_data;
}

void test_CVODE_solver(vector<TestReferenceData> &test_data)
{
    testing::Tester tester = testing::Tester("Test CVODE solver");

    for (auto &data : test_data)
    {
        tester.add_test(data.name, [&]() -> string {
            // Solve the ODE
            OdeFun ode; ode.init(data.cpar);
            OdeSolverCVODE cvode(ode.par->num_species+4);
            OdeSolution solution = cvode.solve(data.t_final, &ode, 60.0, false);
            if (solution.error_ID != ErrorHandler::no_error)
                return ErrorHandler::get_error(solution.error_ID).to_string(true);

            // Calculate molar fractions
            double* x_final = solution.x.back().data();
            vector<double> mol_fractions_final = vector<double>(x_final+3, x_final + ode.par->num_species+3);
            double M = 0.0;
            for (auto &mol_fraction : mol_fractions_final)
                M += mol_fraction;
            for (auto &mol_fraction : mol_fractions_final)
                mol_fraction /= M;

            // Compare the results
            ASSERT_APPROX(solution.x.back().at(2), data.T_final, 1.0e-3);
            for (size_t i = 0; i < data.mol_fractions_final.size(); i++)
            {
                // Check if result is valid number
                if (!std::isfinite(mol_fractions_final[i]))
                {
                    return "mol_fractions_final[" + std::to_string(i) + "] is not a number: " + std::to_string(mol_fractions_final[i]);
                }

                // calculate the difference
                const double abs_diff = std::abs(mol_fractions_final[i] - data.mol_fractions_final[i]);
                const double max = std::max(std::abs(mol_fractions_final[i]), std::abs(data.mol_fractions_final[i]));
                const double rel_diff = abs_diff / max;

                // assemble error message
                std::stringstream ss;
                ss << ode.par->species_names.at(i);
                ss << ": mol_fractions_final[" << i << "] != reference: ";
                ss << std::scientific << std::setprecision(10) << mol_fractions_final[i] << " != " << data.mol_fractions_final[i];
                ss << std::setprecision(3) << " (abs_diff: " << abs_diff << ", rel_diff: " << rel_diff << ")";
                //ss << std::endl << solution;

                // check the difference
                if (max < 1.0e-10)
                {
                    if (abs_diff > 1.0e-12)
                    {
                        return ss.str();
                    }
                }
                else if (max < 1.0e-6)
                {
                    if (rel_diff > 1.0e-2)
                    {
                        return ss.str();
                    }
                }
                else
                {
                    if (rel_diff > 1.0e-3)
                    {
                        return ss.str();
                    }
                }
            }

            return "";
        });
    }

    tester.run_tests();
}



}   // namespace testing

#endif // TEST
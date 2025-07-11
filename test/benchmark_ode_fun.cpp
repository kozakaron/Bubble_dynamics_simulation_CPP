#ifdef BENCHMARK

#include "common.h"
#include "parameters.h"
#include "control_parameters.h"
#include "ode_fun.h"
#include "test_list.h"

// WARNING: THe correctness of the input values are not guaranteed.
void benchmark_ode_fun()
{
    std::cout << colors::bold << "OdeFun() with chemkin_otomo2018 (32 species, 213 reactions)" << colors::reset << std::endl;
    {
    // Set up control parameters
        OdeFun ode = OdeFun();
        const Parameters *par = Parameters::get_parameters(Parameters::mechanism::chemkin_otomo2018);
        (void)par;
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

        (void)ode.init(cpar);

    // Setup inputs
        const double t = 2.5275786980147761e-05;
        std::array<double, 32+4> x = {
            4.35062207e-07, 4.87202821e+00, 4.06585821e+03, 3.04892226e-06,
            3.05074849e-04, 4.22341571e-01, 5.15008532e-07, 2.01019027e-03,
            2.06121865e-06, 8.22242270e-04, 4.39352481e-07, 2.71009207e-01,
            6.41710788e-06, 1.01692955e-05, 4.20519003e-07, 3.61935060e-08,
            3.20915772e-04, 1.37025184e-06, 1.86821312e-07, 1.38126629e-07,
            4.92327628e-07, 2.32080892e-13, 8.92888443e-10, 1.86215771e-08,
            3.03803895e-10, 4.31461142e-15, 1.09808740e-12, 7.88255393e-07,
            7.70016174e-06, 4.84341261e-05, 2.50989113e-05, 2.41829246e-06,
            0.00000000e+00, 0.00000000e+00, 1.40656155e-01, 1.14395306e-06
        };
        std::array<double, 32+4> dxdt;
        const double T = x[2];

    // pressure()
        const double p = 115718.99999999997;
        const double p_dot = 1.4057234685268682e-05;
        std::pair<double, double> result_pressure;
        BENCHMARK_LINE(result_pressure = ode.pressures(t, x[0], x[1], p, p_dot);, 1000000);

    // thermodynamic()
        BENCHMARK_LINE(ode.thermodynamic(T);, 100000);

    // evaporation()
        std::pair<double, double> result_evap;
        const double X_H2O = 3.2315225982636570e-01;
        BENCHMARK_LINE(result_evap = ode.evaporation(p, T, X_H2O);, 1000000);

    // forward_rate()
        double M = 0.8363084533423445;
        std::array<double, 23> M_eff = {4.44204225, 4.44204225, 4.44204225, 1.55048452, 4.77157568,
        3.80254336, 0.83630845, 0.83630845, 3.26860315, 3.26860315,
        0.83630845, 3.26860315, 0.92062359, 0.83630845, 0.83630845,
        0.83630845, 0.83630845, 0.83630845, 0.83630845, 3.90748095,
        0.83630845, 0.83630845, 0.83630845};
        std::copy(M_eff.begin(), M_eff.end(), ode.M_eff);
        
        BENCHMARK_LINE(ode.forward_rate(T, M, p, 1e30);, 100000);

    // backward_rate()
        BENCHMARK_LINE(ode.backward_rate(T, 1e30);, 10000);

    // production_rate()
        BENCHMARK_LINE(ode.production_rate(T, p, M, x.data()+3);, 10000);

    // operator()
        BENCHMARK_LINE(ode(t, x.data(), dxdt.data());, 10000);
    }

    std::cout << colors::bold << "OdeFun() with chemkin_ar_he (12 species, 30 reactions)" << colors::reset << std::endl;
    {
// Set up control parameters
        OdeFun ode = OdeFun();
        ControlParameters cpar;
        const Parameters *par = Parameters::get_parameters(Parameters::mechanism::chemkin_otomo2018);

        cpar.ID = 0;
        cpar.mechanism = Parameters::mechanism::chemkin_ar_he;
        // Initial conditions:
        cpar.R_E = 10e-6;
        cpar.set_species({par->get_species("O2")}, {1.0});
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
        cpar.set_excitation_params({-2.0e5, 30000.0, 1.0});
        cpar.excitation_type = Parameters::excitation::sin_impulse;

        (void)ode.init(cpar);

    // Setup inputs
        const double t = 2.5255799280947053e-05;
        std::array<double, 12+4> x = {
            4.1866326662121154e-07, -1.0806071026205871e-02,  3.7345808358390459e+03,  2.6193400419650303e-05,  2.2272784074252325e-04,
            1.2216575395711056e-03,  6.1721264294022393e-01,  1.9348243023154722e-02,  2.9559164582684627e-01,  0.0000000000000000e+00,
            1.1445537293517848e-02,  5.5500381673767307e-03,  0.0000000000000000e+00,  0.0000000000000000e+00,  7.2195229829469272e-08,
            1.1381442517029756e-06
        };
        std::array<double, 12+4> dxdt;
        const double T = x[2];

    // pressure()
        const double p = 115718.99999999997;
        const double p_dot = 1.4057234685268682e-05;
        std::pair<double, double> result_pressure;
        BENCHMARK_LINE(result_pressure = ode.pressures(t, x[0], x[1], p, p_dot);, 1000000);

    // thermodynamic()
        BENCHMARK_LINE(ode.thermodynamic(T);, 100000);

    // evaporation()
        std::pair<double, double> result_evap;
        const double X_H2O = 3.2315225982636570e-01;
        BENCHMARK_LINE(result_evap = ode.evaporation(p, T, X_H2O);, 1000000);

    // forward_rate()
        double M = 0.8363084533423445;
        std::array<double, 7> M_eff = {
            4.202460954083506 , 4.202460954083506 , 4.202460954083506 , 4.202460954083506 , 4.016609122669752,
            2.0327918805035803, 2.2060452245106026
        };
        std::copy(M_eff.begin(), M_eff.end(), ode.M_eff);
        
        BENCHMARK_LINE(ode.forward_rate(T, M, p, 1e30);, 100000);

    // backward_rate()
        BENCHMARK_LINE(ode.backward_rate(T, 1e30);, 10000);

    // production_rate()
        BENCHMARK_LINE(ode.production_rate(T, p, M, x.data()+3);, 10000);

    // operator()
        BENCHMARK_LINE(ode(t, x.data(), dxdt.data());, 10000);
    }
}

#endif // BENCHMARK
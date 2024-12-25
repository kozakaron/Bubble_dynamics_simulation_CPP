#ifdef BENCHMARK
#ifndef BENCHMARK_ODE_FUN_H
#define BENCHMARK_ODE_FUN_H

#include "ode_fun.h"
#include "test.h"
#include "test_list.h"

// WARNING: THe correctness of the input values are not guaranteed.
void benchmark_ode_fun()
{
// Set up control parameters
    ODE ode = ODE();
    cpar_t cpar;

    cpar.ID = 0;
    cpar.par = Parameters::mechanism::chemkin_otomo2018;
    // Initial conditions:
    cpar.R_E = 10e-6;
    cpar.ratio = 1.0;
    cpar.set_species({2, 31}, {0.75, 0.25});  // TODO: change indexes
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
    cpar.target_specie = 0;//TODO: par::index::NH3;
    // Excitation parameters:
    cpar.set_excitation_params({-2.0e5, 30000.0, 1.0});
    cpar.excitation_type = Parameters::excitation::sin_impulse;

    ode.init(cpar);

// Setup inputs
    const double t = 0.0;
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
    //const double R = x[0];
    //const double R_dot = x[1];
    const double T = x[2];
    std::copy(x.begin(), x.end(), ode.x);

// pressure()
    const double p = 115718.99999999997;
    const double p_dot = 1.4057234685268682e-05;
    std::pair<double, double> result_pressure;
    BENCHMARK_LINE(result_pressure = ode.pressures(t, p, p_dot);, 1000000);

// thermodynamic()
    BENCHMARK_LINE(ode.thermodynamic(T);, 100000);

// evaporation()
    std::pair<double, double> result_evap;
    const double X_H2O = 3.2315225982636570e-01;
    BENCHMARK_LINE(result_evap = ode.evaporation(p, T, X_H2O);, 1000000);

// forward_rate()
    double M = 0.8363084533423445;
    array<double, 23> M_eff = {4.44204225, 4.44204225, 4.44204225, 1.55048452, 4.77157568,
    3.80254336, 0.83630845, 0.83630845, 3.26860315, 3.26860315,
    0.83630845, 3.26860315, 0.92062359, 0.83630845, 0.83630845,
    0.83630845, 0.83630845, 0.83630845, 0.83630845, 3.90748095,
    0.83630845, 0.83630845, 0.83630845};
    std::copy(M_eff.begin(), M_eff.end(), ode.M_eff);
    
    BENCHMARK_LINE(ode.forward_rate(T, M, p);, 100000);

// backward_rate()
    BENCHMARK_LINE(ode.backward_rate(T);, 10000);

// production_rate()
    BENCHMARK_LINE(ode.production_rate(T, p, M);, 10000);
}

#endif // BENCHMARK_ODE_FUN_H
#endif // BENCHMARK
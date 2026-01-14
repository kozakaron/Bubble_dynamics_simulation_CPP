#ifndef ODE_FUN_H
#define ODE_FUN_H

#include "parameters.h"
#include "control_parameters.h"

// Calculate the vapour pressure of water [Pa] as a function of temperature [K]
double vapour_pressure(const double T);

// Calculate the dynamic viscosity of water [Pa*s] as a function of temperature [K] (Pressure dependence is neglected)
double viscosity(const double T);

class OdeFun
{
public:
// Members
    // generic
    const Parameters* par;              // reaction mechanism
    ControlParameters cpar;             // control parameters
    size_t num_species;                 // number of species (to check if init was called properly)
#if defined TEST || defined BENCHMARK
public:
#else
private:
#endif
    // dimensionless form
    double* x_dimensional;              // state vector (x) in SI units (length: Parameters::num_species+4) [mol/m^3]
    // evaporation
    double C_v_inf;                     // molar heat capacity of water at constant volume of ambient temperature [J/mol/K]
    // thermodynamic
    double* C_p;                        // molar heat capacities at constant pressure (length: Parameters::num_species) [J/mol/K]
    double* H;                          // enthalpies (length: Parameters::num_species) [J/mol]
    double* S;                          // entropies (length: Parameters::num_species) [J/mol/K]
    // production rates
    double* M_eff;                      // effective molar masses of third bodies (length: Parameters::num_third_body_reactions)
    double* ln_k_forward;               // logarithm of forward reaction rates (length: Parameters::num_reactions)
    double* ln_k_backward;              // logarithm of backward reaction rates (length: Parameters::num_reactions)
    double* net_rates;                  // net reaction rates (length: Parameters::num_reactions)
    double* omega_dot;                  // production rates (length: Parameters::num_species)

// Methods
    double internal_pressure(
        const double T,
        const double M,
        const double* conc
    );


    double internal_pressure_derivative(
        const double T,
        const double T_dot,
        const double M,
        const double M_dot,
        const double* conc,
        const double* conc_dot
    );


    std::pair<double, double> pressures_excitation(
        const double t,
        const double R,
        const double R_dot,
        const double p,
        const double p_dot
    );


    void thermodynamic(
        const double T
    );


    std::pair<double, double> evaporation(
        const double p,
        const double T,
        const double X_H2O
    );


    void forward_rate(
        const double T,
        const double M,
        const double p
    );


    void backward_rate(
        const double T
    );

    
    void production_rate(
        const double T,
        const double M,
        const double p,
        const double* c
    );

public:
    OdeFun();
    ~OdeFun();
    // Initializes: allocate required memory for current simulation
    is_success init(const ControlParameters& cpar);
    // Calculates the initial condition of the bubble from the control parameters.
    is_success initial_conditions(
        double* x_dimless
    );
    // Call operator. Calculates the right-hand side of the ODE system.
    is_success operator()(
        const double t_dimless,
        const double* x_dimless,
        double* x_dimless_dot
    );

#if defined TEST || defined BENCHMARK
public:
#else
private:
#endif
    void delete_memory();
    is_success check_before_call(const double* x_dimless);
    is_success check_after_call(
        const double t_dimless,
        const double* x_dimless,
        double* x_dimless_dot
    );

};  // class OdeFun

#endif // ODE_FUN_H
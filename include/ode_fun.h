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
    // evaporation
    double C_v_inf;                     // molar heat capacity at constant volume of ambient temperature [erg/mol/K]
    // thermodynamic
    double* C_p;                        // molar heat capacities at constant pressure (length: Parameters::num_species)
    double* H;                          // enthalpies (length: Parameters::num_species)
    double* S;                          // entropies (length: Parameters::num_species)
    double* C_v;                        // molar heat capacities at constant volume (length: Parameters::num_species)
    // production rates
    double* M_eff;                      // effective molar masses of third bodies (length: Parameters::num_third_bodies)
    double* k_forward;                  // forward reaction rates (length: Parameters::num_reactions)
    double* k_backward;                 // backward reaction rates (length: Parameters::num_reactions)
    double* net_rates;                  // net production rates (length: Parameters::num_reactions)
    double* omega_dot;                  // production rates (length: Parameters::num_species)

// Methods
    std::pair<double, double> pressures(
        const double t,
        const double R,
        const double R_dot,
        const double p,
        const double p_dot
    ) ; //noexcept


    void thermodynamic(
        const double T
    ) ; //noexcept


    std::pair<double, double> evaporation(
        const double p,
        const double T,
        const double X_H2O
    ) ; //noexcept


    void forward_rate(
        const double T,
        const double M,
        const double p,
        const double reaction_rate_threshold
    ) ; //noexcept


    void backward_rate(
        const double T,
        const double reaction_rate_threshold
    ) ; //noexcept

    
    void production_rate(
        const double T,
        const double M,
        const double p,
        const double* c
    ) ; //noexcept

public:
    OdeFun();
    ~OdeFun();
    // Initializes: allocate required memory for current simulation
    is_success init(const ControlParameters& cpar);
    // Calculates the initial condition of the bubble from the control parameters.
    is_success initial_conditions(
        double* x
    ) ; //noexcept
    // Call operator. Calculates the right-hand side of the ODE system.
    is_success operator()(
        const double t,
        const double* x,
        double* dxdt
    ) ; //noexcept

#if defined TEST || defined BENCHMARK
public:
#else
private:
#endif
    void delete_memory();
    is_success check_before_call(const double* x);
    is_success check_after_call(
        const double t,
        const double* x,
        double* dxdt
    );

};  // class OdeFun

#endif // ODE_FUN_H
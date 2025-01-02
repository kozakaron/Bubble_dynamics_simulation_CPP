#ifndef ODE_FUN_H
#define ODE_FUN_H

#include "parameters.h"
#include "control_parameters.h"

class ODE
{
private:
    double* memory_block;               // a single allocation for all memory needs
#if defined TEST || defined BENCHMARK
    // WARNING: ODE class manages it's own memory.
    //          These are public for ease of use, however it is not recommended to modify them directly.
    //          Do not modify these pointers directly, use std::copy() instead.
    public:
#else
    private:
#endif
// Members
    // generic
    const Parameters* par;              // reaction mechanism
    cpar_t* cpar;                       // control parameters
    double* x;                          // state vector (length: Parameters::num_species+4)
    double* dxdt;                       // time derivative of state vector (length: Parameters::num_species+4)
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
        const double p
    ) ; //noexcept


    void backward_rate(
        const double T
    ) ; //noexcept

    
    void production_rate(
        const double T,
        const double M,
        const double p
    ) ; //noexcept

public:
    ODE();
    ~ODE();
    void init(const cpar_t& cpar);

private:
    void delete_memory();

};  // class ODE

#endif // ODE_FUN_H
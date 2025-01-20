#ifndef ODE_FUN_H
#define ODE_FUN_H

#include "parameters.h"
#include "control_parameters.h"

class ODE
{
#if defined TEST || defined BENCHMARK
public:
#else
private:
#endif
// Members
    // generic
    const Parameters* par;              // reaction mechanism
    cpar_t* cpar;                       // control parameters
    size_t error_ID;                    // ID of error in ErrorHandler (ErrorHandler::no_error if no error occured)
    size_t num_species;                 // number of species (to check if init was called properly)
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
        const double p
    ) ; //noexcept


    void backward_rate(
        const double T
    ) ; //noexcept

    
    void production_rate(
        const double T,
        const double M,
        const double p,
        const double* c
    ) ; //noexcept

public:
    ODE();
    ~ODE();
    is_success init(const cpar_t& cpar);
    is_success operator()(   // <----- CALL OPERATOR
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
    is_success check_before_call();
    is_success check_after_call(
        const double t,
        const double* x,
        double* dxdt
    );

};  // class ODE

#endif // ODE_FUN_H
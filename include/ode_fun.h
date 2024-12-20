#ifndef ODE_FUN_H
#define ODE_FUN_H
#include <initializer_list>

#include "common.h"
#include "parameters.h"

// ODE function state vector
typedef std::array<double, par::num_species+4> state_t;

/*________________________________control_parameters________________________________*/

typedef class ControlParameters {
public:
// Members
    size_t ID;                  // ID of control parameter
    // Initial conditions:
    double R_E;                 // bubble equilibrium radius [m]
    double ratio;               // R_0/R_E [-]
    par::index* species;        // indexes of species in initial bubble (pointer to array of species index enum)
    double* fractions;          // molar fractions of species in initial bubble (pointer to array of doubles)
    size_t n_species;           // number of species in initial bubble
    // Ambient parameters:
    double P_amb;               // ambient pressure [Pa]
    double T_inf;               // ambient temperature [K]
    // Liquid parameters:
    double alfa_M;              // water accommodation coefficient [-]
    double P_v;                 // vapour pressure [Pa]
    double mu_L;                // dynamic viscosity [Pa*s]
    double rho_L;               // density [kg/m^3]
    double c_L;                 // sound speed [m/s]
    double surfactant;          // surface tension modifier [-]
    // Simulation settings:
    bool enable_heat_transfer;
    bool enable_evaporation;
    bool enable_reactions;
    bool enable_dissipated_energy;
    par::index target_specie;
    // Excitation parameters:
    double* excitation_params;  // parameters for excitation (pointer to array of doubles)
    par::excitation excitation_type; // type of excitation

    
// Methods
    ControlParameters();
    ~ControlParameters();
    // Set species and their fractions like this: set_species({par::index::H2, par::index::N2}, {0.75, 0.25});
    void set_species(const std::initializer_list<par::index> species_list, const std::initializer_list<double> fractions_list);
    // Set excitation parameters like this: set_excitation_params({-2.0e5, 30000.0, 1.0});
    void set_excitation_params(const std::initializer_list<double> params_list);
    // Copy from another ControlParameters object
    void copy(const ControlParameters& cpar);
    // TODO: check()
} cpar_t;

/*________________________________ode_function________________________________*/

class ODE
{
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
    cpar_t* cpar;               // control parameters
    double* x;                  // state vector (length: par::num_species+4)
    // thermodynamic
    double* C_p;                // molar heat capacities at constant pressure (length: par::num_species)
    double* H;                  // enthalpies (length: par::num_species)
    double* S;                  // entropies (length: par::num_species)
    double* C_v;                // molar heat capacities at constant volume (length: par::num_species)
    // production rates
    double* M_eff;              // effective molar masses of third bodies (length: par::num_third_bodies)
    double* k_forward;          // forward reaction rates (length: par::num_reactions)
    double* k_backward;         // backward reaction rates (length: par::num_reactions)

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
        const double T
    ) ; //noexcept


    void forward_rate(
        const double T,
        const double M,
        const double p
    ) ; //noexcept


    void backward_rate(
        const double T
    ) ; //noexcept

public:
    ODE();
    ~ODE();

};  // class ODE

#endif // ODE_FUN_H
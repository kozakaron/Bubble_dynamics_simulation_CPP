#ifndef ODE_FUN_H
#define ODE_FUN_H

#include <tuple>
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
    // Computes internal gas pressure according to ideal gas law or van der Waals state equation. Returns p [Pa].
    double internal_pressure(
        const double T,
        const double M,
        const double* conc
    );


    // Computes the time derivative of internal gas pressure from internal_pressure(). Returns p_dot [Pa/s].
    double internal_pressure_derivative(
        const double T,
        const double T_dot,
        const double M,
        const double M_dot,
        const double* conc,
        const double* conc_dot
    );


    // Computes far field pressure excitation and its time derivative. Returns {P_inf [Pa], P_inf_dot [Pa/s]}.
    std::pair<double, double> excitation_pressures(
        const double t
    );


    // Computes liquid properties at the bubble wall using NASG or Tait EoS. Returns {c_L [m/s], rho_L [kg/m^3], rho_inf [kg/m^3], H [m^2/s^2], T_L [K]}.
    // Returns constant values for Keller-Miksis.
    std::tuple<double, double, double, double, double> liquid_properties(const double p_L, const double P_inf);


    // Computes the second derivative of the bubble radius using Keller-Miksis or Gilmore equation. Returns R_dot_dot [m/s^2].
    double bubble_dynamics(
        const double t,
        const double R,
        const double R_dot,
        const double p,
        const double p_dot,
        const double P_inf,
        const double P_inf_dot
    );


    // Computes thermodynamic properties (heat capacities, enthalpies, entropies) for all species at temperature T using NASA polynomials. Fills internal arrays C_p [J/mol/K], H [J/mol], S [J/mol/K]
    void thermodynamic(
        const double T
    );


    // Computes net molar flux of water into the bubble and associated energy flux. Returns {n_net_dot [mol/m^2/s], evap_energy [W/m^2]}.
    std::pair<double, double> evaporation(
        const double p,
        const double T,
        const double X_H2O
    );


    // Computes forward reaction rates per reaction. Fills internal array ln_k_forward [1/s or m^3/mol/s].
    // supports: Arrhenius, Lindemann, Troe, SRI, PLOG, third-body reactions
    // Can use uni- and bimolecular reaction rate thresholds
    void forward_rate(
        const double T,
        const double M,
        const double p
    );


    // Computes backward reaction rates per reaction. Fills internal array ln_k_backward [1/s or m^3/mol/s].
    void backward_rate(
        const double T
    );

    
    // Computes net production rates per species. Fills internal array omega_dot [mol/m^3/s].
    void production_rate(
        const double T,
        const double M,
        const double p,
        const double* c
    );


    // Call operator. Calculates the right-hand side of the ODE system.
    is_success operator()(
        const double t_dimless,
        const double* x_dimless,
        double* x_dimless_dot
    );


    // Initializes: allocate required memory for current simulation
    is_success init(const ControlParameters& cpar);


    // Calculates the initial condition of the bubble from the control parameters.
    is_success initial_conditions(
        double* x_dimless
    );


    OdeFun();
    ~OdeFun();
    void delete_memory();
    is_success check_before_call(const double* x_dimless);
    is_success check_after_call(
        const double t_dimless,
        const double* x_dimless,
        double* x_dimless_dot
    );

};  // class OdeFun

#endif // ODE_FUN_H
#ifndef CONTROL_PARAMETERS_H
#define CONTROL_PARAMETERS_H
#include <initializer_list>

#include "common.h"
#include "parameters.h"

/*________________________________control_parameters________________________________*/

typedef class ControlParameters {
public:
// Members
    index_t ID;                         // ID of control parameter
    Parameters::mechanism mechanism;    // rection mechanism
    size_t error_ID;                    // ID of error in ErrorHandler (ErrorHandler::no_error if no error occured)
    // Initial conditions:
    double R_E;                         // bubble equilibrium radius [m]
    double ratio;                       // R_0/R_E [-]
    index_t* species;                   // indexes of species in initial bubble (pointer to array of species index enum)
    double* fractions;                  // molar fractions of species in initial bubble (pointer to array of doubles)
    index_t n_species;                  // number of species in initial bubble
    // Ambient parameters:
    double P_amb;                       // ambient pressure [Pa]
    double T_inf;                       // ambient temperature [K]
    // Liquid parameters:
    double alfa_M;                      // water accommodation coefficient [-]
    double P_v;                         // vapour pressure [Pa]
    double mu_L;                        // dynamic viscosity [Pa*s]
    double rho_L;                       // density [kg/m^3]
    double c_L;                         // sound speed [m/s]
    double surfactant;                  // surface tension modifier [-]
    // Simulation settings:
    bool enable_heat_transfer;
    bool enable_evaporation;
    bool enable_reactions;
    bool enable_dissipated_energy;
    index_t target_specie;
    // Excitation parameters:
    double* excitation_params;          // parameters for excitation (pointer to array of doubles)
    Parameters::excitation excitation_type; // type of excitation

    
// Methods
    ControlParameters();
    ~ControlParameters();
    // Set species and their fractions like this: set_species({par->get_species("H2"), par->get_species("N2")}, {0.75, 0.25});
    void set_species(const std::initializer_list<index_t> species_list, const std::initializer_list<double> fractions_list);
    // Set excitation parameters like this: set_excitation_params({-2.0e5, 30000.0, 1.0});
    void set_excitation_params(const std::initializer_list<double> params_list);
    // Copy from another ControlParameters object
    void copy(const ControlParameters& cpar);
} cpar_t;


#endif // CONTROL_PARAMETERS_H
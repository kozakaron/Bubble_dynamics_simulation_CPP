#ifndef CONTROL_PARAMETERS_H
#define CONTROL_PARAMETERS_H
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <initializer_list>

#include "common.h"
#include "parameters.h"

/*________________________________control_parameters________________________________*/

typedef class ControlParameters {
public:
// Constants
    static constexpr size_t max_excitation_params = std::ranges::max(Parameters::excitation_arg_nums);
    static constexpr size_t max_species = 4;
// Members
    size_t ID;                          // ID of control parameter
    Parameters::mechanism mechanism;    // rection mechanism
    size_t error_ID;                    // ID of error in ErrorHandler (ErrorHandler::no_error if no error occured)
    // Initial conditions:
    double R_E;                         // bubble equilibrium radius [m]
    index_t species[max_species];       // indexes of species in initial bubble (pointer to array of species index enum)
    double fractions[max_species];      // molar fractions of species in initial bubble (pointer to array of doubles)
    index_t num_initial_species;        // number of species in initial bubble
    // Ambient parameters:
    double P_amb;                       // ambient pressure [Pa]
    double T_inf;                       // ambient temperature [K]
    // Liquid parameters:
    double alfa_M;                      // water accommodation coefficient [-]
    double P_v;                         // vapour pressure [Pa]
    double mu_L;                        // dynamic viscosity [Pa*s]
    double rho_L;                       // liquid density [kg/m^3]
    double c_L;                         // sound speed [m/s]
    double surfactant;                  // surface tension modifier [-]
    // Simulation settings:
    bool enable_heat_transfer;
    bool enable_evaporation;
    bool enable_reactions;
    bool enable_dissipated_energy;
    index_t target_specie;
    // Excitation parameters:
    double excitation_params[max_excitation_params];    // parameters for excitation (pointer to array of doubles)
    Parameters::excitation excitation_type;             // type of excitation

// Builder struct, defaults
    /* Usage in initialization: cpar_t cpar{cpar_t::Builder{
        .ID = 1,
        .mechanism = Parameters::mechanism::chemkin_ar_he,
        ...
    });
    */
    struct Builder {
        size_t ID                               = 0;
        size_t error_ID                         = ErrorHandler::no_error;
        Parameters::mechanism mechanism         = Parameters::mechanism::chemkin_ar_he;
        double R_E                              = 10.0e-6;
        std::vector<std::string> species        = {"O2"};
        std::vector<double> fractions           = {1.0};
        double P_amb                            = 101325.0;
        double T_inf                            = 293.15;
        double alfa_M                           = 0.35;
        double P_v                              = 2338.1;
        double mu_L                             = 0.001;
        double rho_L                            = 998.2;
        double c_L                              = 1483.0;
        double surfactant                       = 1.0;
        bool enable_heat_transfer               = true;
        bool enable_evaporation                 = true;
        bool enable_reactions                   = true;
        bool enable_dissipated_energy           = true;
        std::string target_specie               = "H2";
        std::vector<double> excitation_params   = {-2.0e5, 30000.0, 1.0};
        Parameters::excitation excitation_type  = Parameters::excitation::sin_impulse;
    };
    
// Methods
    ControlParameters();
    ControlParameters(const Builder& builder);
    ~ControlParameters();
    // Set species and their fractions like this: set_species({"H2", "N2"}, {0.75, 0.25}); or set_species({par->get_species("O2")}, {1.0});
    void set_species(const std::vector<std::string> species_list, const std::vector<double> fractions_list);
    void set_species(const std::vector<index_t>& species_list, const std::vector<double>& fractions_list);
    void set_species(const std::initializer_list<std::string>& species_list, const std::initializer_list<double>& fractions_list);
    void set_species(const std::initializer_list<index_t>& species_list, const std::initializer_list<double>& fractions_list);
    // Set excitation parameters like this: set_excitation_params({-2.0e5, 30000.0, 1.0});
    void set_excitation_params(const std::vector<double>& params_list);
    void set_excitation_params(const std::initializer_list<double>& params_list);
    // Convert to string. with_code: ready to copy to constructor; one_line: without newlines
    std::string to_string(const bool with_code=false) const;
    // Ostream overload
    friend std::ostream& operator<<(std::ostream& os, const ControlParameters& cpar);
private:
    void init(const Builder& builder);
} cpar_t;


#endif // CONTROL_PARAMETERS_H
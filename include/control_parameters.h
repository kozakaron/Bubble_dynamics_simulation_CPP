#ifndef CONTROL_PARAMETERS_H
#define CONTROL_PARAMETERS_H
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <initializer_list>

#include "nlohmann/json_fwd.hpp"
#include "common.h"
#include "parameters.h"


class ControlParameters {
public:
// Constants
    static constexpr size_t max_excitation_params = std::ranges::max(Parameters::excitation_arg_nums);
    static constexpr size_t max_species = 5;
    static constexpr char csv_header[] = "ID,mechanism,R_E,species,fractions,P_amb,T_inf,alfa_M,P_v,mu_L,rho_L,c_L,surfactant,enable_heat_transfer,enable_evaporation,enable_reactions,enable_dissipated_energy,target_specie,excitation_params,excitation_type";
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
    // Reference values for dimensionless parameters: A_star = A / A_ref
    static constexpr double t_ref = 1e-6;               // reference time [s]
    static constexpr double t_ref_inv = 1e6;            // inverse reference time [1/s]
    double R_ref;                                       // reference radius [m]
    double T_ref;                                       // reference temperature [K]
    static constexpr double E_diss_ref = 1e-5;          // reference dissipated energy [J] (instead of R_ref^2 * t_ref_inv^2)
    static constexpr double epsilon = 1e-30;            // c_star = ln(c + epsilon)


// Builder struct, defaults
    /* Usage in initialization: ControlParameters cpar{ControlParameters::Builder{
        .ID = 1,
        .mechanism = Parameters::mechanism::chemkin_ar_he,
        ...
    }};
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
    // Constructor with Builder. See Builder struct for more details.
    ControlParameters(const Builder& builder);
    // Constructor from JSON (ordered_json). Input should contains that same keys as in Builder struct, but mechanism and excitation_type are strings.
    // Example: {"ID": 1, "mechanism": "chemkin_ar_he", "species": ["O2"], "fractions": [1.0], ...}
    ControlParameters(const nlohmann::ordered_json& j);
    // Constructor from JSON file (using 'cpar' key)
    ControlParameters(const std::string& json_path);
    ~ControlParameters();
    // Set species and their fractions like this: set_species({"H2", "N2"}, {0.75, 0.25}); or set_species({par->get_species("O2")}, {1.0});
    void set_species(const std::vector<std::string> species_list, const std::vector<double> fractions_list);
    void set_species(const std::vector<index_t>& species_list, const std::vector<double>& fractions_list);
    void set_species(const std::initializer_list<std::string>& species_list, const std::initializer_list<double>& fractions_list);
    void set_species(const std::initializer_list<index_t>& species_list, const std::initializer_list<double>& fractions_list);
    // Set excitation parameters like this: set_excitation_params({-2.0e5, 30000.0, 1.0});
    void set_excitation_params(const std::vector<double>& params_list);
    void set_excitation_params(const std::initializer_list<double>& params_list);
    // Convert to CSV string
    std::string to_csv() const;
    // Convert to string. with_code: ready to copy to constructor; one_line: without newlines
    std::string to_string(const bool with_code=false) const;
    // Convert to JSON
    nlohmann::ordered_json to_json() const;
    // Ostream overload
    friend std::ostream& operator<<(std::ostream& os, const ControlParameters& cpar);
    // Transform from x(t) to x_star(t_star)
    void nondimensionalize(double &t, double* x) const;
    // Transform from x_star(t_star) to x(t)
    void dimensionalize(double &t, double* x) const;
private:
    void init(const Builder& builder);
};


#endif // CONTROL_PARAMETERS_H
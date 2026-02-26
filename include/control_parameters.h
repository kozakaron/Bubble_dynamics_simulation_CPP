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
    static constexpr size_t max_species = 10;
    static constexpr char csv_header[] = "ID,mechanism,R_E,ratio,species,fractions,P_amb,T_inf,alpha_M,P_v,mu_L,rho_L,c_L,surfactant,enable_heat_transfer,enable_evaporation,enable_reactions,enable_dissipated_energy,enable_van_der_waals,enable_rate_thresholding,target_specie,excitation_type,excitation_params,excitation_cycles,ramp_up_cycles";
// Members
    size_t ID;                          // ID of control parameter
    const Parameters* par;              // reaction mechanism (pointer to Parameters instance)
    size_t error_ID;                    // ID of error in ErrorHandler (ErrorHandler::no_error if no error occured)
    // Initial conditions:
    double R_E;                         // bubble equilibrium radius [m]
    double ratio;                       // R_0/R_E for unforced oscillations [-]
    index_t species[max_species];       // indexes of species in initial bubble (array of species index enum)
    double fractions[max_species];      // molar fractions of species in initial bubble (array of doubles)
    index_t num_initial_species;        // number of species in initial bubble
    // Ambient parameters:
    double P_amb;                       // ambient pressure [Pa]
    double T_inf;                       // ambient temperature [K]
    // Liquid parameters:
    double alpha_M;                     // water accommodation coefficient [-]
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
    bool enable_van_der_waals;
    bool enable_rate_thresholding;
    index_t target_specie;
    
    size_t R_dot_from_file;
    size_t rows;
    size_t cols;
	std::string file_name;
	std::vector<double> importdata;
    
    // Excitation parameters:
    Parameters::excitation excitation_type;             // type of excitation
    double excitation_params[max_excitation_params];    // parameters for excitation (pointer to array of doubles)
    double excitation_cycles;                           // number of excitation cycles to use (according to freq/freq1 in excitation_params) [-]
    double ramp_up_cycles;                              // number of cycles until the excitation reaches full amplitude (0<=ramp_up_cycles<=excitation_cycles/2) [-]

    // Reference values for dimensionless parameters: A_dimless = A / A_ref
    static constexpr double t_ref = 1e-9;               // reference time [s]
    static constexpr double t_ref_inv = 1e9;            // inverse reference time [1/s]
    double n_ref;                                       // reference concentration [mol/m3]
    double R_ref;                                       // reference radius [m]
    double T_ref;                                       // reference temperature [K]
    double E_diss_ref;                                  // reference dissipated energy [J] (instead of R_ref^2 * t_ref_inv^2)
	
	//Reference values for Gilmore:
	/*double Gamma_L;
	double c_L_ref;    
	double kappa;
	double r_hc;
	double B_L;
	double b_L;
	double cV_L;
	double rho_L_ref;
	double p_L_ref;*/

// Builder struct, defaults
    /* Usage in initialization: ControlParameters cpar{ControlParameters::Builder{
        .ID = 1,
        .mechanism = "chemkin_elte2016_hydrogen",
        ...
    }};
    */
    struct Builder {
        size_t ID                               = 0;
        size_t error_ID                         = ErrorHandler::no_error;
        std::string mechanism                   = "chemkin_elte2016_hydrogen";
        double R_E                              = 10.0e-6;
        double ratio                            = 1.0;
        std::vector<std::string> species        = {"O2"};
        std::vector<double> fractions           = {1.0};
        double P_amb                            = 101325.0;
        double T_inf                            = 293.15;
        double alpha_M                          = 0.35;
        double P_v                              = 2338.1;
        double mu_L                             = 0.001;
        double rho_L                            = 998.2;
        double c_L                              = 1483.0;
        double surfactant                       = 1.0;
		/*double Gamma_L							= 1.19;
		double c_L_ref							= 1496.0;
		double kappa							= 1.4;
		double r_hc								= 3.537973231411522e-06;
		double B_L								= 616640000.0;
		double b_L								= 6.72e-4;
		double cV_L								= 3657.19;
		double rho_L_ref						= 998.2;
		double p_L_ref							= 1.0e5;		*/
        bool enable_heat_transfer               = true;
        bool enable_evaporation                 = true;
        bool enable_reactions                   = true;
        bool enable_dissipated_energy           = true;
        bool enable_van_der_waals               = true;
        bool enable_rate_thresholding           = true;

		size_t R_dot_from_file					= 0;
		size_t rows 							= 0;
		size_t cols 							= 0;
		std::string file_name						= "";
		std::vector<double> importdata;
		
        std::string target_specie               = "H2";
        Parameters::excitation excitation_type  = Parameters::excitation::sinusoid;
        std::vector<double> excitation_params   = {-2.0e5, 30000.0};
        double excitation_cycles                = 1.0;
        double ramp_up_cycles                   = 0.0;
    };
    
// Methods
    ControlParameters();
    // Constructor with Builder. See Builder struct for more details.
    ControlParameters(const Builder& builder);
    // Constructor from JSON (ordered_json). Input should contains that same keys as in Builder struct, but mechanism and excitation_type are strings.
    // Example: {"ID": 1, "mechanism": "chemkin_elte2016_hydrogen", "species": ["O2"], "fractions": [1.0], ...}
    ControlParameters(const nlohmann::ordered_json& j);
    // Constructor from JSON file (using 'cpar' key)
    ControlParameters(const std::string& json_path);
    ~ControlParameters();
    // Set reaction mechanism by name
    void set_mechanism(const std::string& mechanism_name);
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
    // Transform from x(t) to x_dimless(t_dimless) in place
    void nondimensionalize(double &t, double* x) const;
    // Transform from x_dimless(t_dimless) to x(t) in place
    void dimensionalize(double &t, double* x) const;
    // Transform x_dot to x_dimless_dot in place
    void nondimensionalize_dot(double* x_dot, const double* x) const;
private:
    void init(const Builder& builder);
};


#endif // CONTROL_PARAMETERS_H
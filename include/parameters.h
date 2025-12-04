#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <array>
#include <vector>
#include <string>
#include <unordered_map>
#include <cstdint>
#include <memory>

#include "nlohmann/json_fwd.hpp"

typedef uint16_t index_t;   // Type for indexes (like size_t)
typedef int8_t stoich_t;    // Type for stoichiometric coefficients (type of constants in reactions: 2*H2 + O2 -> 2*H2O)

class Parameters
{
private:
    static std::vector<std::unique_ptr<Parameters>> _mechanisms;     // Vector for already loaded mechanisms
    std::unordered_map<std::string, index_t> _species;               // Names of species (num_species)
public:
// EXCITATION AND REACTION TYPES

    enum reac_type: index_t {lindemann, troe, sri};
    enum excitation: index_t {no_excitation=0, sinusoid=1, two_sinusoids=2, square=3};
    static constexpr std::array<index_t, 4> excitation_arg_nums = {
        0, // no_excitation
        2, // sinusoid
        5, // two_sinusoids
        3  // square
    };
    static constexpr std::array<const char*, 4> excitation_names = {
        "no_excitation",
        "sinusoid",
        "two_sinusoids",
        "square"
    };
    static constexpr std::array<const char*, 4> excitation_arg_names = {
        "",                                            // no_excitation
        "p_A freq",                                    // sinusoid
        "p_A1 p_A2 freq1 freq2 phase_shift",           // two_sinusoids
        "p_A freq harmonics"                           // square
    };
    static constexpr std::array<const char*, 4> excitation_arg_units = {
        "",                                            // no_excitation
        "Pa Hz",                                       // sinusoid
        "Pa Pa Hz Hz rad",                             // two_sinusoids
        "Pa Hz -"                                      // square
    };

// PHYSICAL CONSTANTS

    static constexpr double c_L           = 1483.0;             // Liquid sound speed at 30 °C [m/s]
    static constexpr double rho_L         = 998.2;              // Liquid density [kg/m^3]
    static constexpr double sigma         = 0.07197;            // Surface tension [N/m]
    static constexpr double mu_L          = 0.001;              // Dynamic viscosity at 30 °C and 1 atm [Pa*s]
    static constexpr double P_v           = 2338.1;             // Saturated vapour pressure at 30 °C [Pa]
    static constexpr double alpha_M       = 0.35;               // Water accommodation coefficient [-]
    static constexpr double k_B           = 1.380649e-23;       // Boltzmann constant [J/K]
    static constexpr double R_g           = 8.31446;            // Universal gas constant [J/mol/K]
    static constexpr double R_kmol        = 8314.46;            // Universal gas constant [J/kmol/K]
    static constexpr double R_erg         = 83144600.0;         // Universal gas constant [erg/mol/K]
    static constexpr double R_cal         = 1.987204;           // Universal gas constant [cal/mol/K]
    static constexpr double N_A           = 6.02214e+23;        // Avogadro's number [-]
    static constexpr double ln_N_A        = 54.754899816742695;  // Natural logarithm of Avogadro's number
    static constexpr double h             = 6.62607015e-34;     // Planck constant [m^2*kg/s]
    static constexpr double R_v           = 461.521126;         // Specific gas constant of water [J/kg/K]
    static constexpr double erg2J         = 1e-07;              // Conversion factor from erg to J
    static constexpr double cal2J         = 4.184;              // Conversion factor from cal to J
    static constexpr double atm2Pa        = 101325.0;           // Conversion factor from atm to Pa
    static constexpr double bar2Pa        = 100000.0;           // Conversion factor from bar to Pa
    static constexpr double absolute_zero = 273.15;             // Zero °C in Kelvin

// MECHANISM DEPENDENT PARAMETERS

// Common mechanism data
    const std::string mechanism_name;               // Name of the mechanism
// Species
    const index_t num_elements;                     // Number of elements
    const index_t num_species;                      // Number of species
    const index_t index_of_water;                   // Index of water in arrays, INVALID if H2O is not in mechanism
    const index_t invalid_index;                    // Invalid index
    std::vector<std::string> species_names;         // Names of species (num_species)
    index_t get_species(std::string name) const;    // Get index of species by name
    const double *W;                                // Molar masses [kg/mol] (num_species)
    const double *thermal_conductivities;           // Thermal conductivities [W/m/K] (num_species)
    const double *sqrt_van_der_waals_a;             // square root of van der Waals a parameters [m^3*Pa^0.5/mol] (num_species)
    const double *van_der_waals_b;                  // van der Waals b parameters [m^3/mol] (num_species)
// NASA polynomials
    static constexpr index_t NASA_order=5;          // Degree of NASA polynomials
    const double *temp_ranges;                      // Temperature ranges for NASA polynomials [K] (num_species, 3)
    const double *a_low;                            // Low  NASA coefficients (num_species, NASA_order+2)
    const double *a_high;                           // High NASA coefficients (num_species, NASA_order+2)
    const double *interval_values;                  // Values at T_low and T_high (num_species, 3): {C_p_high, H_high, S_high}
    const double *interval_derivatives;             // Derivatives at T_low and T_high (num_species, 3): {dC_p_high/dT, dH_high/dT, dS_high/dT}
// Reactions constants
    const index_t num_reactions;                    // Number of reactions
    const double *arrhenius_parameters;             // Arrhenius parameters (num_reactions, 3): {ln_A, b, E_over_R}
    const index_t *reaction_order;                  // reaction_order (num_reactions)
// Reaction matrixes
    const index_t num_max_species_per_reaction;     // Maximum number of species participating in a reaction
    const index_t *nu_indexes;                      // Indexes of species participating in reactions (num_reactions, num_max_species_per_reaction)
    const stoich_t *nu_forward;                     // Stoichiometric coefficients left hand side (num_reactions, num_max_species_per_reaction)
    const stoich_t *nu_backward;                    // Stoichiometric coefficients right hand side (num_reactions, num_max_species_per_reaction)
    const stoich_t *nu;                             // Stoichiometric coefficients nu = nu_backward - nu_forward (num_reactions, num_max_species_per_reaction)
    const stoich_t *sum_nu;                         // Sum of stoichiometric coefficients for each reaction (num_reactions)
// Third body reactions
    const index_t num_third_body_reactions;         // Number of third body reactions
    const index_t *third_body_reaction_indexes;     // Indexes of third body species (num_third_body_reactions)
    const bool *is_falloff_reaction;                // Pressure dependent flag (num_third_body_reactions)
    const double *alpha;                            // Third body efficiency factors (num_third_body_reactions, num_species)
// Irreversible reactions
    const index_t num_irreversible_reactions;       // Number of irreversible reactions
    const index_t *irreversible_reaction_indexes;   // Indexes of irreversible reactions (num_irreversible_reactions)
// Pressure-dependent fall-off reactions
    const index_t num_falloff_reactions;            // Number of pressure dependent reactions
    const index_t num_lindemann_reactions;          // Number of Lindemann reactions
    const index_t num_troe_reactions;               // Number of Troe reactions
    const index_t num_sri_reactions;                // Number of SRI reactions
    const index_t *falloff_reaction_indexes;        // Indexes of pressure dependent reactions (num_falloff_reactions)
    const reac_type *falloff_reaction_types;        // Types of pressure dependent reactions (num_falloff_reactions)
    const index_t *is_third_body_indexes;           // Indexes of third body species for pressure dependent reactions (num_falloff_reactions)
    const double *falloff_parameters;               // Fall-off parameters (num_pressure_dependent, 3): {ln_A0, b_0, E0_over_R}
    const double *troe_parameters;                  // Troe parameters (num_troe_reactions, 4): {alpha, T***, T**, T*}
    const double *sri_parameters;                   // SRI parameters (num_sri_reactions, 5): {a, b, c, d, e}
// Pressure-dependent PLOG reactions
    const index_t num_plog_reactions;               // Number of PLOG reactions
    const index_t num_plog_levels;                  // Number of PLOG levels
    const index_t *plog_reaction_indexes;           // Indexes of PLOG reactions (num_plog_reactions)
    const index_t *plog_seperators;                 // Seperators of PLOG reactions (num_plog_reactions+1)
    const double *plog_parameters;                  // PLOG parameters (num_plog_levels, 4): {P_j, ln_Aj, b_j, Ej_over_R}

// CONSTRUCTORS
    Parameters(const nlohmann::json& j);
    Parameters(const std::string& json_path);
    ~Parameters();

// GETTERS
    static const Parameters *get_parameters(const std::string& mech_name);
    static Parameters::excitation string_to_excitation(std::string excitation_str);
};

#endif // PARAMETERS_H
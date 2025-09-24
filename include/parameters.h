#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <array>
#include <vector>
#include <string>
#include <unordered_map>
#include <cstdint>

typedef uint16_t index_t;   // Type for indexes (like size_t)
typedef int8_t stoich_t;    // Type for stoichiometric coefficients (type of constants in reactions: 2*H2 + O2 -> 2*H2O)

class Parameters
{
private:
    static const Parameters chemkin_ar_he_params;
    static const Parameters chemkin_kaust2023_n2_params;
    static const Parameters chemkin_kaust2023_n2_without_o_params;
    static const Parameters chemkin_otomo2018_without_o_params;
    static const Parameters chemkin_otomo2018_params;
    static const Parameters chemkin_elte2016_ethanol_params;
    static const Parameters chemkin_elte2016_syngas_params;
    static const Parameters chemkin_elte2017_methanol_params;

    std::unordered_map<std::string, index_t> _elements;       // Names of elements (num_elements)
    std::unordered_map<std::string, index_t> _species;        // Names of species (num_species)
public:
// MECHANISM, EXCITATION AND REACTION TYPES

    enum mechanism: index_t {
        chemkin_ar_he,
        chemkin_kaust2023_n2,
        chemkin_kaust2023_n2_without_o,
        chemkin_otomo2018_without_o,
        chemkin_otomo2018,
        chemkin_elte2016_ethanol,
        chemkin_elte2016_syngas,
        chemkin_elte2017_methanol
    };
    enum reac_type: index_t {lindemann_reac, troe_reac, sri_reac};
    enum excitation: index_t {no_excitation=0, two_sinusoids=1, sin_impulse=2, sin_impulse_logf=3};
    static constexpr std::array<index_t, 4> excitation_arg_nums = {
        0, // no_excitation
        6, // two_sinusoids
        3, // sin_impulse
        3  // sin_impulse_logf
    };
    static constexpr std::array<const char*, 4> excitation_names = {
        "no_excitation",
        "two_sinusoids",
        "sin_impulse",
        "sin_impulse_logf"
    };
    static constexpr std::array<const char*, 4> excitation_arg_names = {
        "",                                            // no_excitation
        "p_A1 p_A2 freq1 freq2 theta_phase n_cycles",  // two_sinusoids
        "p_A freq n_cycles",                           // sin_impulse
        "p_A log_f n_cycles"                           // sin_impulse_logf
    };
    static constexpr std::array<const char*, 4> excitation_arg_units = {
        "",                                            // no_excitation
        "Pa Pa Hz Hz rad -",                           // two_sinusoids
        "Pa Hz -",                                     // sin_impulse
        "Pa - -"                                       // sin_impulse_logf
    };
    static constexpr std::array<const char*, 8> mechanism_names = {
        "chemkin_ar_he",
        "chemkin_kaust2023_n2",
        "chemkin_kaust2023_n2_without_o",
        "chemkin_otomo2018_without_o",
        "chemkin_otomo2018",
        "chemkin_elte2016_ethanol",
        "chemkin_elte2016_syngas",
        "chemkin_elte2017_methanol"
    };

// PHYSICAL CONSTANTS

    static constexpr double c_L           = 1483.0;             // Liquid sound speed at 30 °C [m/s]
    static constexpr double rho_L         = 998.2;              // Liquid density [kg/m^3]
    static constexpr double sigma         = 0.07197;            // Surface tension [N/m]
    static constexpr double mu_L          = 0.001;              // Dynamic viscosity at 30 °C and 1 atm [Pa*s]
    static constexpr double P_v           = 2338.1;             // Saturated vapour pressure at 30 °C [Pa]
    static constexpr double alfa_M        = 0.35;               // Water accommodation coefficient [-]
    static constexpr double k_B           = 1.380649e-23;       // Boltzmann constant [J/K]
    static constexpr double R_g           = 8.31446;            // Universal gas constant [J/mol/K]
    static constexpr double R_erg         = 83144600.0;         // Universal gas constant [erg/mol/K]
    static constexpr double R_cal         = 1.987204;           // Universal gas constant [cal/mol/K]
    static constexpr double N_A           = 6.02214e+23;        // Avogadro's number [-]
    static constexpr double log_N_A       = 54.754899816742695; // Natural logarithm of Avogadro's number
    static constexpr double h             = 6.62607015e-34;     // Planck constant [m^2*kg/s]
    static constexpr double R_v           = 461.521126;         // Specific gas constant of water [J/kg/K]
    static constexpr double erg2J         = 1e-07;              // Conversion factor from erg to J
    static constexpr double cal2J         = 4.184;              // Conversion factor from cal to J
    static constexpr double atm2Pa        = 101325.0;           // Conversion factor from atm to Pa
    static constexpr double bar2Pa        = 100000.0;           // Conversion factor from bar to Pa
    static constexpr double absolute_zero = 273.15;             // Zero °C in Kelvin

// MECHANISM DEPENDENT PARAMETERS

// Common mechanism data
    const std::string model;                        // Name of the mechanism
    const std::string input_file;                   // Original input (.inp) file's name
// Species
    const index_t num_elements;                     // Number of elements
    const index_t num_species;                      // Number of species
    const index_t index_of_water;                   // Index of water in arrays, INVALID if H2O is not in mechanism
    const index_t invalid_index;                    // Invalid index
    std::vector<std::string> elements_names;        // Names of elements (num_elements)
    std::vector<std::string> species_names;         // Names of species (num_species)
    index_t get_element(std::string name) const;    // Get index of element by name
    index_t get_species(std::string name) const;    // Get index of species by name
    const double *W;                                // Molar masses [g/mol] (num_species)
    const double *lambdas;                          // Thermal conductivities [W/m/K] (num_species)
// NASA polynomials
    const index_t NASA_order;                       // Degree of NASA polynomials
    const double *temp_range;                       // Temperature ranges for NASA polynomials [K] (num_species, 3)
    const double *a_low;                            // Low  NASA coefficients (num_species, NASA_order+2)
    const double *a_high;                           // High NASA coefficients (num_species, NASA_order+2)
    const double *interval_values;                  // Values at T_low and T_high (num_species, 3): {C_p_high, H_high, S_high}
    const double *interval_derivatives;             // Derivatives at T_low and T_high (num_species, 3): {dC_p_high/dT, dH_high/dT, dS_high/dT}
// Reactions constants
    const index_t num_reactions;                    // Number of reactions
    const double *b;                                // Temperature exponents [-] (num_reactions)
    const double *ln_A;                             // Logarithm of pre-exponential factors [ln(cm^3/mol/s) v ln(1/s)] (num_reactions)
    const double *E_over_R;                         // Activation energies / universal gas constant [K] (num_reactions)
    const double *N_A_pow_reaction_order;           // N_A^reaction_order (num_reactions)
// Reaction matrixes
    const index_t num_max_specie_per_reaction;      // Maximum number of species participating in a reaction
    const index_t *nu_indexes;                      // Indexes of species participating in reactions (num_reactions, num_max_specie_per_reaction)
    const stoich_t *nu_forward;                     // Stoichiometric coefficients left hand side (num_reactions, num_max_specie_per_reaction)
    const stoich_t *nu_backward;                    // Stoichiometric coefficients right hand side (num_reactions, num_max_specie_per_reaction)
    const stoich_t *nu;                             // Stoichiometric coefficients nu = nu_backward - nu_forward (num_reactions, num_max_specie_per_reaction)
    const stoich_t *sum_nu;                         // Sum of stoichiometric coefficients for each reaction (num_reactions)
// Third body reactions
    const index_t num_third_bodies;                 // Number of third body reactions
    const index_t *third_body_indexes;              // Indexes of third body species (num_third_bodies)
    const bool *is_pressure_dependent;              // Pressure dependent flag (num_third_bodies)
    const double *alfa;                             // Third body efficiency factors (num_third_bodies, num_species)
// Irreversible reactions
    const index_t num_irreversible;                 // Number of irreversible reactions
    const index_t *irreversible_indexes;            // Indexes of irreversible reactions (num_irreversible)
// Pressure dependent reactions
    const index_t num_pressure_dependent;           // Number of pressure dependent reactions
    const index_t num_lindemann;                    // Number of Lindemann reactions
    const index_t num_troe;                         // Number of Troe reactions
    const index_t num_sri;                          // Number of SRI reactions
    const index_t *pressure_dependent_indexes;      // Indexes of pressure dependent reactions (num_pressure_dependent)
    const reac_type *pressure_dependent_reac_types; // Types of pressure dependent reactions (num_pressure_dependent)
    const index_t *is_third_body_indexes;           // Indexes of third body species for pressure dependent reactions (num_pressure_dependent)
    const double *reac_const;                       // Fall-off parameters (num_pressure_dependent, 3): {ln_A0, b_0, E0_over_R}
    const double *troe;                             // Troe parameters (num_troe, 4): {alfa, T***, T**, T*}
    const double *sri;                              // SRI parameters (num_sri, 5): {a, b, c, d, e}
// PLOG reactions
    const index_t num_plog;                         // Number of PLOG reactions
    const index_t num_plog_levels;                  // Number of PLOG levels
    const index_t *plog_indexes;                    // Indexes of PLOG reactions (num_plog)
    const index_t *plog_seperators;                 // Seperators of PLOG reactions (num_plog+1)
    const double *plog;                             // PLOG parameters (num_plog_levels, 4): {P_j, ln_Aj, b_j, Ej_over_R}

// CONSTRUCTORS
    template <typename Parameters_struct>
    Parameters(Parameters_struct dummy);
    ~Parameters();

// GETTERS
    static const Parameters *get_parameters(const Parameters::mechanism mech);
    static Parameters::mechanism string_to_mechanism(std::string mechanism_str);
    static Parameters::excitation string_to_excitation(std::string excitation_str);
};

#endif // PARAMETERS_H
#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <array>
#include <string>

using std::array;
using std::string;

namespace par{

enum reac_type {LINDEMANN, TROE, SRI};
enum excitation {no_excitation=0, two_sinusoids=1, sin_impulse=2, sin_impulse_logf=3};
constexpr array<size_t, 4> excitation_arg_nums = {
    0, // no_excitation
    5, // two_sinusoids
    3, // sin_impulse
    3  // sin_impulse_logf
};
const array<string, 4> excitation_arg_names = {
    "",                                     // no_excitation
    "p_A1 p_A2 freq1 freq2 theta_pahse",    // two_sinusoids
    "p_A freq n",                           // sin_impulse
    "p_A log_f n"                           // sin_impulse_logf
};
const array<string, 4> excitation_arg_units = {
    "",                                     // no_excitation
    "Pa Pa Hz Hz rad",                      // two_sinusoids
    "Pa Hz -",                              // sin_impulse
    "Pa - -"                                // sin_impulse_logf
};


/*________________________________Physical constants________________________________*/

constexpr double c_L           = 1483.0;             // Liquid sound speed at 30 °C [m/s]
constexpr double rho_L         = 998.2;              // Liquid density [kg/m^3]
constexpr double sigma         = 0.07197;            // Surface tension [N/m]
constexpr double mu_L          = 0.001;              // Dynamic viscosity at 30 °C and 1 atm [Pa*s]
constexpr double P_v           = 2338.1;             // Saturated vapour pressure at 30 °C [Pa]
constexpr double alfa_M        = 0.35;               // Water accommodation coefficient [-]
constexpr double R_g           = 8.31446;            // Universal gas constant [J/mol/K]
constexpr double R_erg         = 83144600.0;         // Universal gas constant [erg/mol/K]
constexpr double R_cal         = 1.987204;           // Universal gas constant [cal/mol/K]
constexpr double N_A           = 6.02214e+23;        // Avogadro's number [-]
constexpr double h             = 6.62607015e-34;     // Planck constant [m^2*kg/s]
constexpr double R_v           = 461.521126;         // Specific gas constant of water [J/kg/K]
constexpr double erg2J         = 1e-07;              // Conversion factor from erg to J
constexpr double cal2J         = 4.184;              // Conversion factor from cal to J
constexpr double atm2Pa        = 101325.0;           // Conversion factor from atm to Pa
constexpr double bar2Pa        = 100000.0;           // Conversion factor from bar to Pa
constexpr double absolute_zero = 273.15;             // Zero °C in Kelvin


}   // namespace par

#ifdef CHEMKIN_AR_HE
#include "chemkin_ar_he.h"
#elif defined CHEMKIN_OTOMO2018_WITHOUT_O
#include "chemkin_otomo2018_without_o.h"
#elif defined CHEMKIN_OTOMO2018
#include "chemkin_otomo2018.h"
#else
#error No mechanism defined
#endif

#endif // PARAMETERS_H
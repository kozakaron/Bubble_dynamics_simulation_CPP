#include <algorithm>
#include <sstream>
#include <iomanip>
#include <numbers>
#include <cmath>
#include <tuple>

#include "common.h"
#include "ode_fun.h"


// Computes saturated vapour pressure of water at temperature T [K] using the Arden-Buck equation.
double vapour_pressure(const double T)
{
    double T_C = T - 273.15; // [°C]
    return 611.21 * std::exp( (18.678 - T_C / 234.5) * (T_C / (257.14 + T_C)) ); // [Pa]
}



// Computes dynamic viscosity of water [Pa*s] at temperature T [K] using 4 parameter exponential formula.
double viscosity(const double T)
{
    return 1.856e-14 * std::exp(4209.0/T + 0.04527*T - 3.376e-5*T*T); // [Pa*s]
}


OdeFun::OdeFun():
    par(nullptr),
    cpar(),
    num_species(0),
    x_dimensional(nullptr),
    C_p(nullptr),
    H(nullptr),
    S(nullptr),
    M_eff(nullptr),
    ln_k_forward(nullptr),
    ln_k_backward(nullptr),
    net_rates(nullptr),
    omega_dot(nullptr)
{ }


OdeFun::~OdeFun()
{
    this->delete_memory();
}


void OdeFun::delete_memory()
{
    if (this->x_dimensional != nullptr)  delete[] this->x_dimensional;
    if (this->C_p != nullptr)            delete[] this->C_p;
    if (this->H != nullptr)              delete[] this->H;
    if (this->S != nullptr)              delete[] this->S;
    if (this->M_eff != nullptr)          delete[] this->M_eff;
    if (this->ln_k_forward != nullptr)   delete[] this->ln_k_forward;
    if (this->ln_k_backward != nullptr)  delete[] this->ln_k_backward;
    if (this->net_rates != nullptr)      delete[] this->net_rates;
    if (this->omega_dot != nullptr)      delete[] this->omega_dot;

    this->x_dimensional = nullptr;
    this->C_p = nullptr;
    this->H = nullptr;
    this->S = nullptr;
    this->M_eff = nullptr;
    this->ln_k_forward = nullptr;
    this->ln_k_backward = nullptr;
    this->net_rates = nullptr;
    this->omega_dot = nullptr;
}

is_success OdeFun::check_before_call(const double* x_dimless)
{
    // Check if initialization was correct
    if (this->cpar.error_ID != ErrorHandler::no_error)
    {
        return false;
    }
    else if (this->par == nullptr || this->omega_dot == nullptr)
    {
        this->cpar.error_ID = LOG_ERROR(Error::severity::error, Error::type::odefun, "Arrays/par were nullptr. Forgot to call OdeFun::init()?", this->cpar.ID);
        return false;
    }
    else if (this->num_species != this->par->num_species)
    {
        std::string message = "Invalid array lengths. Expected num_species=" + std::to_string(this->par->num_species) +\
                              ", arrays are initialized with size " + std::to_string(this->num_species) + " instead. Forgot to call OdeFun::init()?";
        this->cpar.error_ID = LOG_ERROR(Error::severity::error, Error::type::odefun, message, this->cpar.ID);
        return false;
    }
    
    // Check if R and T are valid
    if (x_dimless == nullptr)
    {
        this->cpar.error_ID = LOG_ERROR(Error::severity::error, Error::type::odefun, "State vector x_dimless is nullptr", this->cpar.ID);
        return false;
    }
    else if(!std::isfinite(x_dimless[0]) || !std::isfinite(x_dimless[2]))
    {
        std::string message = "Non finite R or T: R=" + std::to_string(x_dimless[0]) + ";   T=" + std::to_string(x_dimless[2]);
        LOG_ERROR(Error::severity::warning, Error::type::odefun, message, this->cpar.ID);  // recoverable error
        return false;
    }
    else if(x_dimless[0] < 0 || x_dimless[2] < 0)
    {
        std::string message = "Negative R or T: R=" + std::to_string(x_dimless[0]) + ";   T=" + std::to_string(x_dimless[2]);
        LOG_ERROR(Error::severity::warning, Error::type::odefun, message, this->cpar.ID);  // recoverable error
        return false;
    }
    
    return true;
}


is_success OdeFun::check_after_call(
    const double t_dimless,
    const double* x_dimless,
    double* x_dimless_dot
) {
    for (index_t k = 0; k < this->num_species+4; ++k)
    {
        if (!std::isfinite(x_dimless_dot[k]))
        {
            std::stringstream ss;
            if (std::isinf(x_dimless_dot[k]))
            {
                ss << "x_dimless_dot[" << k << "] is infinite";
            }
            else if (std::isnan(x_dimless_dot[k]))
            {
                ss << "x_dimless_dot[" << k << "] is NaN";
            } else
            {
                ss << "x_dimless_dot[" << k << "] is not finite";
            }
            ss << ". t_dimless = " << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10) << t_dimless;
            ss << ";    x_dimless = " << to_string((double*)x_dimless, this->num_species+4);
            ss << ";    x_dimless_dot = " << to_string((double*)x_dimless_dot, this->num_species+4);
            this->cpar.error_ID = LOG_ERROR(Error::severity::error, Error::type::odefun, ss.str(), this->cpar.ID);
            return false;
        }
    }
    return true;
}


is_success OdeFun::init(const ControlParameters& cpar)
{
    this->cpar = cpar;
    const Parameters *old_par = this->par;
    this->par = cpar.par;
    if (this->par == nullptr) return false;
    this->num_species = this->par->num_species;
    if (this->cpar.enable_evaporation && this->par->index_of_water == this->par->invalid_index)
    {
        LOG_ERROR(Error::severity::warning, Error::type::odefun, "Water is not in the mechanism, but evaporation is enabled.", this->cpar.ID);
    }
    
    if (old_par != this->par || this->omega_dot == nullptr)
    {
        this->delete_memory();
        this->x_dimensional = new double[par->num_species+4];
        this->C_p           = new double[par->num_species];
        this->H             = new double[par->num_species];
        this->S             = new double[par->num_species];
        this->M_eff         = new double[par->num_third_body_reactions];
        this->ln_k_forward  = new double[par->num_reactions];
        this->ln_k_backward = new double[par->num_reactions];
        this->net_rates     = new double[par->num_reactions];
        this->omega_dot     = new double[par->num_species];
    }
    if(this->cpar.error_ID != ErrorHandler::no_error)
    {
        return false;
    }
    
    // Evaporation calculations
    if (this->par->index_of_water == this->par->invalid_index)
    {
        this->cpar.enable_evaporation = false;
    }
    if (!this->cpar.enable_evaporation)
    {
        this->C_v_inf = 0.0;
        return true;
    }

    // calculate C_v_inf for water
    const double& T_mid = par->temp_ranges[3*par->index_of_water+1];
    const double& T_high = par->temp_ranges[3*par->index_of_water+2];
    
    if (this->cpar.T_inf < T_high)
    {
    // get NASA coefficients: a={a1, a2, a3, a4, a5, a6, a7}
        const double *a;
        if (this->cpar.T_inf < T_mid)
            a = &(par->a_low[par->index_of_water*(par->NASA_order+2)]);
        else
            a = &(par->a_high[par->index_of_water*(par->NASA_order+2)]);

    // calculate polynomial
        const double T1 = this->cpar.T_inf;
        const double T2 = T1 * T1;
        const double T3 = T2 * T1;
        const double T4 = T2 * T2;

        this->C_v_inf = par->R_g * (a[0] + a[1]*T1 + a[2]*T2 + a[3]*T3 + a[4]*T4 - 1.0);
    }
    else
    {
    // linear extrapolation
        const double *interval_values      = &(par->interval_values[par->index_of_water*3]);       // {C_p_high, H_high, S_high}
        const double *interval_derivatives = &(par->interval_derivatives[par->index_of_water*3]);  // {dC_p_high/dT, dH_high/dT, dS_high/dT}
        this->C_v_inf = interval_values[0] + interval_derivatives[0] * (this->cpar.T_inf - T_high) - par->R_g;
    }

    // Check excitation parameters (skip if using excitation interpolator)
    if (this->cpar.excitation_interpolator.x_data.empty() && this->cpar.radius_interpolator.x_data.empty())
    {
        if (this->cpar.excitation_cycles <= 0.0)
        {
            this->cpar.error_ID = LOG_ERROR(Error::severity::error, Error::type::odefun, "Number of excitation cycles must be positive, instead it is " + std::to_string(this->cpar.excitation_cycles), this->cpar.ID);
            return false;
        }
        if (this->cpar.ramp_up_cycles < 0.0 || this->cpar.ramp_up_cycles > this->cpar.excitation_cycles / 2.0)
        {
            this->cpar.error_ID = LOG_ERROR(Error::severity::error, Error::type::odefun, "Number of ramp-up cycles must be in the range [0, " + std::to_string(this->cpar.excitation_cycles / 2.0) + "], instead it is " + std::to_string(this->cpar.ramp_up_cycles), this->cpar.ID);
            return false;
        }
    }
    switch (this->cpar.excitation_type)
    {
    case Parameters::excitation::no_excitation:
        break;
    case Parameters::excitation::sinusoid:
        if (this->cpar.excitation_params[1] <= 0.0)
        {
            this->cpar.error_ID = LOG_ERROR(Error::severity::error, Error::type::odefun, "freq must be positive, instead it is " + std::to_string(this->cpar.excitation_params[1]), this->cpar.ID);
            return false;
        }
        break;
    case Parameters::excitation::two_sinusoids:
        if (this->cpar.excitation_params[2] <= 0.0)
        {
            this->cpar.error_ID = LOG_ERROR(Error::severity::error, Error::type::odefun, "freq1 must be positive, instead it is " + std::to_string(this->cpar.excitation_params[2]), this->cpar.ID);
            return false;
        }
        if (this->cpar.excitation_params[3] <= 0.0)
        {
            this->cpar.error_ID = LOG_ERROR(Error::severity::error, Error::type::odefun, "freq2 must be positive, instead it is " + std::to_string(this->cpar.excitation_params[3]), this->cpar.ID);
            return false;
        }
        break;
    case Parameters::excitation::square:
        if (this->cpar.excitation_params[1] <= 0.0)
        {
            this->cpar.error_ID = LOG_ERROR(Error::severity::error, Error::type::odefun, "Frequency must be positive, instead it is " + std::to_string(this->cpar.excitation_params[1]), this->cpar.ID);
            return false;
        }
        if (this->cpar.excitation_params[2] < 1.0 || std::abs(std::round(this->cpar.excitation_params[2]) - this->cpar.excitation_params[2]) > 1e-12)
        {
            this->cpar.error_ID = LOG_ERROR(Error::severity::error, Error::type::odefun, "Number of harmonics must be a positive integer, instead it is " + std::to_string(this->cpar.excitation_params[2]), this->cpar.ID);
            return false;
        }
        break;
    default:
        this->cpar.error_ID = LOG_ERROR(Error::severity::error, Error::type::odefun, "Invalid excitation type enum value: " + std::to_string(this->cpar.excitation_type), this->cpar.ID);
        return false;
    }

    return true;
}


is_success OdeFun::initial_conditions(
    double* x_dimless
)
{
// Van der Waals mixture
    double a_tot = 0.0;
    double b_tot = 0.0;
    if (cpar.enable_van_der_waals)
    {
        for (index_t k = 0; k < cpar.max_species; ++k)
        {
            index_t index = cpar.species[k];
            double X_k = cpar.fractions[k];
            if (index == par->invalid_index) continue;
            a_tot += X_k * par->sqrt_van_der_waals_a[index];
            b_tot += X_k * par->van_der_waals_b[index];
            // fraction of water is ignored
        }
        a_tot = a_tot * a_tot;
    }   

// Equilibrium state pressure
    const double p_E = cpar.P_amb + 2.0 * cpar.surfactant * par->sigma / cpar.R_E;   // [Pa]
    const double V_E = 4.0 / 3.0 * std::numbers::pi * cpar.R_E * cpar.R_E * cpar.R_E;    // [m^3]
    double p_gas = cpar.enable_evaporation ? p_E - cpar.P_v : p_E;   // [Pa]

// Compute molar concentration from pressure
    double c_gas = p_gas / (par->R_g * cpar.T_inf); // [mol/m^3]
    if (cpar.enable_van_der_waals)
    {
        // Equation: f(c_gas) = R_g*T_inf*c_gas / (1 - b_tot * c_gas) - a_tot * c_gas^2 - p_gas = 0
        // Initial guess: ideal gas law
        // Newton-Raphson method
        const int max_iter = 20;
        const double tol = 1e-8;

        for (int i = 0; i < max_iter; ++i)
        {
            double den = 1.0 - c_gas * b_tot;
            if (den < 1e-9) den = 1e-9; // avoid singularity

            double p_calc = (c_gas * par->R_g * cpar.T_inf) / den - a_tot * c_gas * c_gas;
            double df_dc = (par->R_g * cpar.T_inf) / (den * den) - 2.0 * a_tot * c_gas;
            double diff = p_calc - p_gas;
            if (std::abs(diff) < tol) break;

            c_gas = c_gas - diff / df_dc;
        }
    }

// Equilibrium state moles
    const double n_gas = c_gas * V_E;    // [mol]
    double n_H2O = cpar.enable_evaporation ? cpar.P_v * V_E / (par->R_g * cpar.T_inf) : 0.0;   // [mol]
    double c_H2O = n_H2O / V_E;   // [mol/m^3]

// Isotermic expansion
    const double R_0 = cpar.ratio * cpar.R_E;   // [m]
    const double V_0 = 4.0 / 3.0 * std::numbers::pi * R_0 * R_0 * R_0;   // [m^3]
    n_H2O = cpar.enable_evaporation ? cpar.P_v * V_0 / (par->R_g * cpar.T_inf) : 0.0; // [mol]
    c_H2O = n_H2O / V_0;   // [mol/m^3]
    c_gas = n_gas / V_0;   // [mol/m^3]
    p_gas = c_gas * par->R_g * cpar.T_inf; // [Pa]
    const double P_amb_min = (cpar.enable_evaporation ? cpar.P_v : 0.0) + p_gas - 2.0 * cpar.surfactant * par->sigma / R_0; // [Pa]
    if (P_amb_min < (cpar.enable_evaporation ? cpar.P_v : 0.0))
    {
        LOG_ERROR(Error::severity::warning, Error::type::odefun, "The pressure during the expansion is lower, than the saturated water pressure.", this->cpar.ID);
    }
    
// Initial conditions
    double* x_dimensional = x_dimless;
    x_dimensional[0] = R_0;         // R_0 [m]
    x_dimensional[1] = 0.0;         // R_dot_0 [m/s]
    x_dimensional[2] = cpar.T_inf;  // T_0 [K]
    for (index_t k = 0; k < par->num_species; ++k)
    {
        x_dimensional[3+k] = 0.0;   // c_k_0 [mol/m^3]
    }
    x_dimensional[3+par->num_species] = 0.0;   // dissipated energy [J]
    if (cpar.enable_evaporation && par->index_of_water != par->invalid_index)
    {
        x_dimensional[3+par->index_of_water] = c_H2O; // c_H2O_0 [mol/m^3]
    }
    for (index_t k = 0; k < cpar.num_initial_species; ++k)
    {
        index_t index = cpar.species[k];
        x_dimensional[3+index] = cpar.fractions[k] * c_gas;   // c_k_0 [mol/m^3]
    }

// Dimensionless form
    double dummy = 0.0;
    cpar.nondimensionalize(dummy, x_dimensional);
    // x_dimensional --> x_dimless

// Errors
    if (p_gas < 0.0)
    {
        this->cpar.error_ID = LOG_ERROR(Error::severity::error, Error::type::odefun, "Negative gas pressure: " + std::to_string(p_gas), this->cpar.ID);
        return false;
    }
    return true;
}


double OdeFun::internal_pressure(
    const double T,
    const double M,
    const double* conc
)
{
    if (cpar.enable_van_der_waals)
    {
        double A_tot = 0.0;     // a_tot * M^2
        double B_tot = 0.0;     // b_tot * M
        for (index_t k = 0; k < par->num_species; ++k)
        {
            A_tot += conc[k] * par->sqrt_van_der_waals_a[k];
            B_tot += conc[k] * par->van_der_waals_b[k];
        }
        A_tot = A_tot * A_tot;

        return (M * par->R_g * T) / (1.0 - B_tot) - A_tot;
    }
    else
    {
        return M * par->R_g * T;
    }
}


double OdeFun::internal_pressure_derivative(
    const double T,
    const double T_dot,
    const double M,
    const double M_dot,
    const double* conc,
    const double* conc_dot
)
{
    if (cpar.enable_van_der_waals)
    {
        double A_tot     = 0.0;    // a_tot * M^2
        double A_tot_dot = 0.0;
        double B_tot     = 0.0;    // b_tot * M
        double B_tot_dot = 0.0;
        for (index_t k = 0; k < par->num_species; ++k)
        {
            A_tot     += conc[k]     * par->sqrt_van_der_waals_a[k];
            A_tot_dot += conc_dot[k] * par->sqrt_van_der_waals_a[k];
            B_tot     += conc[k]     * par->van_der_waals_b[k];
            B_tot_dot += conc_dot[k] * par->van_der_waals_b[k];
        }
        A_tot_dot = 2.0 * A_tot * A_tot_dot;
        A_tot = A_tot * A_tot;

        const double nom     = par->R_g * T * M;
        const double nom_dot = par->R_g * (T_dot * M + T * M_dot);
        const double den     = 1.0 - B_tot;
        const double den_dot = -B_tot_dot;

        return (nom_dot * den - nom * den_dot) / (den * den) - A_tot_dot;
    }
    else
    {
        return M * par->R_g * T_dot + M_dot * par->R_g * T;
    }
}


// Helper function for Tukey window calculation used in OdeFun::pressures()
inline std::pair<double, double> _tukey_window(
    const double t,
    const double freq,
    const double excitation_cycles,
    const double ramp_up_cycles
)
{
    const double t_end = excitation_cycles / freq;
    const double t_ramp = ramp_up_cycles / freq;

    double p, p_dot;
    if (t < 0.0 || t > t_end)
    {
        p = 0.0;
        p_dot = 0.0;
    }
    else if (t < t_ramp)
    {
        p = 0.5 * (1.0 - cos(std::numbers::pi * t / t_ramp));
        p_dot = 0.5 * (std::numbers::pi / t_ramp) * sin(std::numbers::pi * t / t_ramp);
    }
    else if (t > t_end - t_ramp)
    {
        p = 0.5 * (1.0 - cos(std::numbers::pi * (t_end - t) / t_ramp));
        p_dot = 0.5 * (std::numbers::pi / t_ramp) * sin(std::numbers::pi * (t_end - t) / t_ramp);
    }
    else
    {
        p = 1.0;
        p_dot = 0.0;
    }
    return std::make_pair(p, p_dot);
}


// Helper function for sinusoidal pressure calculation used in OdeFun::pressures()
inline std::pair<double, double> _sin(
    const double t,
    const double p_A,
    const double freq,
    const double phase_shift
)
{
    const double insin = 2.0 * std::numbers::pi * freq;
    const double p = p_A * sin(insin * t + phase_shift);
    const double p_dot = p_A * insin * cos(insin * t + phase_shift);
    return std::make_pair(p, p_dot);
}


std::pair<double, double> OdeFun::excitation_pressures(
    const double t
)
{
    double p_inf, p_inf_dot;
       
    switch (cpar.excitation_type)
    {
    case Parameters::excitation::sinusoid:
        {
            const double& p_A = cpar.excitation_params[0];
            const double& freq = cpar.excitation_params[1];

            const auto [w, w_dot] = _tukey_window(t, freq, cpar.excitation_cycles, cpar.ramp_up_cycles);
            const auto [p1, p1_dot] = _sin(t, p_A, freq, 0.0);
            p_inf = cpar.P_amb + w * p1;
            p_inf_dot = w_dot * p1 + w * p1_dot;
            break;
        }
    case Parameters::excitation::two_sinusoids:
        {
            const double& p_A1 = cpar.excitation_params[0];
            const double& p_A2 = cpar.excitation_params[1];
            const double& freq1 = cpar.excitation_params[2];
            const double& freq2 = cpar.excitation_params[3];
            const double& theta_phase = cpar.excitation_params[4];

            const auto [w, w_dot] = _tukey_window(t, freq1, cpar.excitation_cycles, cpar.ramp_up_cycles);
            const auto [p1, p1_dot] = _sin(t, p_A1, freq1, 0.0);
            const auto [p2, p2_dot] = _sin(t, p_A2, freq2, theta_phase);
            p_inf = cpar.P_amb + w * (p1 + p2);
            p_inf_dot = w_dot * (p1 + p2) + w * (p1_dot + p2_dot);
            break;
        }
    case Parameters::excitation::square:
        {
            const double& p_A = cpar.excitation_params[0];
            const double& freq = cpar.excitation_params[1];
            const size_t harmonics = round(cpar.excitation_params[2]);

            const auto [w, w_dot] = _tukey_window(t, freq, cpar.excitation_cycles, cpar.ramp_up_cycles);
            double p_sum = 0.0;
            double p_sum_dot = 0.0;
            for (size_t n = 1; n <= 2*harmonics; n += 2)
            {
                const double p_An = p_A * (4.0 / (std::numbers::pi * n));
                const double freqn = n * freq;
                const auto [p_n, p_n_dot] = _sin(t, p_An, freqn, 0.0);
                p_sum += p_n;
                p_sum_dot += p_n_dot;
            }
            p_inf = cpar.P_amb + w * p_sum;
            p_inf_dot = w_dot * p_sum + w * p_sum_dot;
            break;
        }
    default: // no_excitation
        {
            p_inf = cpar.P_amb;
            p_inf_dot = 0.0;
            break;
        }
    }

    // Use excitation interpolator if available
    if (!cpar.excitation_interpolator.x_data.empty())
    {
        auto [p_interp, p_dot_interp, p_dot_dot_interp] = cpar.excitation_interpolator.interpolate(t);
        p_inf = cpar.P_amb + p_interp;
        p_inf_dot = p_dot_interp;
    }

    return std::make_pair(p_inf, p_inf_dot);
}


std::tuple<double, double, double, double, double> OdeFun::liquid_properties(const double p_L, const double P_inf)
{
    if (cpar.bubble_dynamics_type == Parameters::bubble_dynamics::keller_miksis)
    {
        return {cpar.c_L, cpar.rho_L, cpar.rho_L, 0.0, cpar.T_inf};
    } else if (cpar.bubble_dynamics_type == Parameters::bubble_dynamics::gilmore_nasg)
    {
        // Noble-Abel stiffened-gas (NASG) EoS
        // Validation of parameter count happens in ControlParameters::init()
        const double Gamma_L   = cpar.liquid_eos_params[0];
        const double B_L       = cpar.liquid_eos_params[1];
        const double b_L       = cpar.liquid_eos_params[2];
        const double p_L_ref   = cpar.liquid_eos_params[3];
        const double rho_L_ref = cpar.liquid_eos_params[4];

        const double K_L       = rho_L_ref / (std::pow(p_L_ref + B_L, 1.0 / Gamma_L) * (1.0 - b_L * rho_L_ref));
        const double rho_L     = K_L * std::pow((p_L   + B_L), 1.0 / Gamma_L) / (1.0 + b_L * K_L * std::pow((p_L   + B_L), 1.0 / Gamma_L));
        const double rho_inf   = K_L * std::pow((P_inf + B_L), 1.0 / Gamma_L) / (1.0 + b_L * K_L * std::pow((P_inf + B_L), 1.0 / Gamma_L));
        const double c_L       = std::sqrt(Gamma_L * (p_L + B_L) / (rho_L - b_L * rho_L * rho_L));
        const double h_L       = Gamma_L / (Gamma_L - 1.0) * (p_L   + B_L) / rho_L   - Gamma_L / (Gamma_L - 1.0) * b_L * (p_L   + B_L) + b_L * p_L;
        const double h_inf     = Gamma_L / (Gamma_L - 1.0) * (P_inf + B_L) / rho_inf - Gamma_L / (Gamma_L - 1.0) * b_L * (P_inf + B_L) + b_L * P_inf;
        const double H         = h_L - h_inf;
        const double T_L       = cpar.T_inf * std::pow((p_L + B_L) / (cpar.P_amb + B_L), (Gamma_L - 1.0) / Gamma_L);

        return {c_L, rho_L, rho_inf, H, T_L};
    } else {
        // Tait EoS (default for Gilmore)
        // Validation of parameter count happens in ControlParameters::init()
        const double Gamma_L   = cpar.liquid_eos_params[0];
        const double B_L       = cpar.liquid_eos_params[1];
        const double p_L_ref   = cpar.liquid_eos_params[2];
        const double rho_L_ref = cpar.liquid_eos_params[3];

        const double rho_L     = rho_L_ref * std::pow((p_L   + B_L) / (p_L_ref + B_L), 1.0 / Gamma_L);
        const double rho_inf   = rho_L_ref * std::pow((P_inf + B_L) / (p_L_ref + B_L), 1.0 / Gamma_L);
        const double c_L       = std::sqrt(Gamma_L * (p_L + B_L) / rho_L);
        const double h_L       = Gamma_L / (Gamma_L - 1.0) * (p_L   + B_L) / rho_L;
        const double h_inf     = Gamma_L / (Gamma_L - 1.0) * (P_inf + B_L) / rho_inf;
        const double H         = h_L - h_inf;
        const double Gamma_G   = 0.1 + 1.4 * (1.0 - rho_L_ref / rho_L);
        const double rho_ratio = rho_L_ref / rho_L;
        const double T_L       = cpar.T_inf * std::exp(Gamma_G * (1.0 - rho_ratio));

        return {c_L, rho_L, rho_inf, H, T_L};
    }
}


double OdeFun::bubble_dynamics(
    const double t,
    const double R,
    const double R_dot,
    const double p,
    const double p_dot,
    const double P_inf,
    const double P_inf_dot
)
{
    // Radius interpolator takes precedence over equation-based dynamics
    if (!cpar.radius_interpolator.x_data.empty())
    {
        auto [R_interp, R_dot_interp, R_dot_dot_interp] = cpar.radius_interpolator.interpolate(t);
        const double damping_coefficient = 1.0e8;    // adjustable [1/s]
        return R_dot_dot_interp - damping_coefficient * (R_dot - R_dot_interp);
    }

    if (cpar.bubble_dynamics_type == Parameters::bubble_dynamics::keller_miksis)
    {
        // Keller-Miksis equation
        const double p_L     = p - (2.0 * cpar.surfactant * par->sigma + 4.0 * cpar.mu_L * R_dot) / R;
        const double p_L_dot = p_dot + (2.0 * cpar.surfactant * par->sigma * R_dot + 4.0 * cpar.mu_L * R_dot * R_dot) / (R * R); // - 4.0 * cpar.mu_L * R_dot_dot / R;
        const double delta = (p_L - P_inf) / cpar.rho_L;
        const double delta_dot = (p_L_dot - P_inf_dot) / cpar.rho_L;

        const double nom   = (1.0 + R_dot / cpar.c_L) * delta + R / cpar.c_L * delta_dot - (1.5 - 0.5 * R_dot / cpar.c_L) * R_dot * R_dot;
        const double denom = (1.0 - R_dot / cpar.c_L) * R + 4.0 * cpar.mu_L / (cpar.c_L * cpar.rho_L);
        return nom / denom;
    }
    else    // Gilmore equation (NASG or Tait EoS) - https://doi.org/10.1016/j.ultsonch.2020.105307
    {
        const double p_L     = p - (2.0 * cpar.surfactant * par->sigma + 4.0 * cpar.mu_L * R_dot) / R;
        const double p_L_dot = p_dot + (2.0 * cpar.surfactant * par->sigma * R_dot + 4.0 * cpar.mu_L * R_dot * R_dot) / (R * R); // - 4.0 * cpar.mu_L * R_dot_dot / R;
        auto [c_L, rho_L, rho_inf, H, T_L_unused] = liquid_properties(p_L, P_inf);
        const double H_dot = p_L_dot / rho_L - P_inf_dot / rho_inf; // - 4.0 * cpar.mu_L * R_dot_dot / (rho_L * R);

        const double nom = ((1.0 + R_dot / c_L) * H - 1.5 * (1.0 - R_dot / (3.0 * c_L)) * R_dot * R_dot) / ((1.0 - R_dot / c_L) * R) + H_dot / c_L;
        const double den = 1.0 + 4.0 * cpar.mu_L / (rho_L * R * c_L);
        return nom / den;
    }
}


void OdeFun::thermodynamic(
    const double T
)
{
// calculate powers for polynomials
    const double T1 = T;
    const double T2 = T1 * T1;
    const double T3 = T2 * T1;
    const double T4 = T2 * T2;
    const double T5 = T3 * T2;

    for (index_t k = 0; k < par->num_species; ++k)
    {
        const double& T_low = par->temp_ranges[3*k];
        const double& T_mid = par->temp_ranges[3*k+1];
        const double& T_high = par->temp_ranges[3*k+2];
        (void)T_low;

        if (T < T_high)
        {
        // get NASA coefficients: a={a1, a2, a3, a4, a5, a6, a7}
            const double *a;
            if (T < T_mid)
            {
                a = &(par->a_low[k*(par->NASA_order+2)]);
            }
            else
            {
                a = &(par->a_high[k*(par->NASA_order+2)]);
            }
        
            // Molar heat capacities at constant pressure (isobaric) [J/mol/K]
            this->C_p[k] = par->R_g * (a[0] + a[1]*T1 + a[2]*T2 + a[3]*T3 + a[4]*T4);
            // Enthalpies [J/mol]
            this->H[k] = par->R_g * (a[0]*T1 + a[1]*T2/2.0 + a[2]*T3/3.0 + a[3]*T4/4.0 + a[4]*T5/5.0 + a[5]);
            // Entropies [J/mol/K]
            this->S[k] = par->R_g * (a[0]*std::log(T) + a[1]*T1 + a[2]*T2/2.0 + a[3]*T3/3.0 + a[4]*T4/4.0 + a[6]);
        }
        else
        {
        // linear extrapolation
            const double *interval_values      = &(par->interval_values[k*3]);       // {C_p_high, H_high, S_high}
            const double *interval_derivatives = &(par->interval_derivatives[k*3]);  // {dC_p_high/dT, dH_high/dT, dS_high/dT}
            this->C_p[k] = interval_values[0] + interval_derivatives[0] * (T - T_high);
            this->H[k]   = interval_values[1] + interval_derivatives[1] * (T - T_high);
            this->S[k]   = interval_values[2] + interval_derivatives[2] * (T - T_high);
        }
    }
}


std::pair<double, double> OdeFun::evaporation(
    const double p,
    const double T,
    const double X_H2O
)
{
// condensation and evaporation
    const double p_H2O = p * X_H2O;
    const double n_eva_dot = cpar.alpha_M * cpar.P_v / (par->W[par->index_of_water] * std::sqrt(2.0 * std::numbers::pi * par->R_v * cpar.T_inf));
    const double n_con_dot = cpar.alpha_M * p_H2O    / (par->W[par->index_of_water] * std::sqrt(2.0 * std::numbers::pi * par->R_v * T));
    const double n_net_dot = n_eva_dot - n_con_dot;
// Evaporation energy [J/mol]
    const double& C_v = this->C_p[par->index_of_water] - par->R_g; // Molar heat capacity of water at constant volume (isochoric) [J/mol/K]
    double e_eva = this->C_v_inf * cpar.T_inf;
    double e_con = C_v * T;
    double evap_energy = n_eva_dot * e_eva - n_con_dot * e_con;    // [W/m^2]

    return std::make_pair(n_net_dot, evap_energy);
}


// Natural logarithm of the dimensionless collision integral Ω^(2,2)(T*).
// https://doi.org/10.1063/1.1678363
inline double ln_Omega_22(const double T_star, const double ln_T_star)
{
    constexpr double A    = 1.16145, B = 0.14874;
    constexpr double C    = 0.52487, D = 0.77320;
    constexpr double E    = 2.16178, F = 2.43787;
    constexpr double ln_A = 0.1496692245;
    if (T_star > 10.0)
    {
        return ln_A - B * ln_T_star;
    }
    return std::log(
        A * std::pow(T_star, -B) + C * std::exp(-D * T_star) + E * std::exp(-F * T_star)
    );
}

void OdeFun::forward_rate(
    const double T,
    const double M,
    const double p
)
{
// Arrhenius reactions
    const double ln_T = std::log(T);
    for(index_t index = 0; index < par->num_reactions; ++index)
    {
        const double& ln_A     = par->arrhenius_parameters[index * 3 + 0];    // Logarithm of pre-exponential factors [ln(m^3n/mol^n/s)], where n is reaction_order-1 [-]
        const double& b        = par->arrhenius_parameters[index * 3 + 1];    // Temperature exponents [-]
        const double& E_over_R = par->arrhenius_parameters[index * 3 + 2];    // Activation energies / universal gas constant [K]

        this->ln_k_forward[index] = ln_A + b * ln_T - E_over_R / T;
    }
    
// Pressure dependent reactions
    index_t troe_index=0, sri_index=0;
    for(index_t j = 0; j < par->num_falloff_reactions; ++j)
    {
    // Third body reactions
        double M_eff_loc = M;
        index_t third_body_index = par->is_third_body_indexes[j];
        if (third_body_index != par->invalid_index)
        {
            M_eff_loc = this->M_eff[third_body_index];
        }
        
    // Pressure dependent formalisms
        index_t index = par->falloff_reaction_indexes[j];
        const double *falloff_parameters = &(par->falloff_parameters[j*3]);
        const double &ln_A0     = falloff_parameters[0];
        const double &b_0       = falloff_parameters[1];
        const double &E0_over_R = falloff_parameters[2];

        const double ln_k_inf = this->ln_k_forward[index];
        const double ln_k_0 = ln_A0 + b_0 * ln_T - E0_over_R / T;
        const double ln_Pr = ln_k_0 + std::log(M_eff_loc) - ln_k_inf;
        const double log10_Pr = ln_Pr / std::numbers::ln10;

        double ln_F = 0.0;  // F = 1.0

        switch (par->falloff_reaction_types[j])
        {
        case Parameters::reac_type::lindemann:
            {
                // F = 1.0;
                break;
            }
        case Parameters::reac_type::troe:
            {
                const double *troe_parameters = &(par->troe_parameters[troe_index*4]);
                const double &alpha  = troe_parameters[0];
                const double &T_xxx = troe_parameters[1];  // T***
                const double &T_x   = troe_parameters[2];  // T*
                const double &T_xx  = troe_parameters[3];  // T**

                constexpr double small_num = 1.0e-30;
                constexpr double large_num = 1.0e30;
                const double exp_negT_over_Txxx = (T_xxx <= small_num) ? 0.0 : std::exp(-T / T_xxx);
                const double exp_negT_over_Tx   = (T_x   >= large_num) ? 1.0 : std::exp(-T / T_x);
                const double exp_negTxx_over_T  = (T_xx  >= large_num) ? 0.0 : std::exp(-T_xx / T);
                const double F_cent = (1.0 - alpha) * exp_negT_over_Txxx + alpha * exp_negT_over_Tx + exp_negTxx_over_T;
                const double log10_F_cent = std::log10(F_cent);
                
                const double c = -0.4 - 0.67 * log10_F_cent;
                const double n = 0.75 - 1.27 * log10_F_cent;
                constexpr double d = 0.14;
                const double e = log10_Pr + c;
                const double log10_F = log10_F_cent / (1.0 + std::pow(e / (n - d * e), 2));
                
                ln_F = log10_F * std::numbers::ln10;
                ++troe_index;
                break;
            }
        case Parameters::reac_type::sri:   // TODO: We don't have any SRI reactions in the current mechanisms
            {
                const double *sri_parameters = &(par->sri_parameters[sri_index*5]);
                const double &a = sri_parameters[0];
                const double &b = sri_parameters[1];
                const double &c = sri_parameters[2];
                const double &d = sri_parameters[3];
                const double &e = sri_parameters[4];

                const double X = 1.0 / (1.0 + log10_Pr * log10_Pr);
                const double F = d * std::pow(a * std::exp(-b / T) + std::exp(-T / c), X) * std::pow(T, e);
                ln_F = std::log(F);
                ++sri_index;
                break;
            }
        default:
            break;
        }

        this->ln_k_forward[index] = ln_k_inf + ln_Pr - std::log1p(std::exp(ln_Pr)) + ln_F;  // std::log(1.0 + P_r) = std::log1p(P_r)
    } // pressure dependent reactions end
    
// PLOG reactions
    for(index_t j = 0; j < par->num_plog_reactions; ++j)
    {
        index_t index = par->plog_reaction_indexes[j];
        // determne indexes of the lower and upper pressures
        index_t lower = par->plog_seperators[j];
        for (index_t k = par->plog_seperators[j] + 1; k < par->plog_seperators[j+1] - 1; ++k)
        {
            const double &P_j = par->plog_parameters[k*4+0];
            if (P_j < p)
                lower = k;
            else
                break;
        }
        index_t upper = lower + 1;
        const double &P_lower        = par->plog_parameters[lower*4+0];
        const double &P_upper        = par->plog_parameters[upper*4+0];
        const double &ln_A_lower     = par->plog_parameters[lower*4+1];
        const double &ln_A_upper     = par->plog_parameters[upper*4+1];
        const double &b_lower        = par->plog_parameters[lower*4+2];
        const double &b_upper        = par->plog_parameters[upper*4+2];
        const double &E_lower_over_R = par->plog_parameters[lower*4+3];
        const double &E_upper_over_R = par->plog_parameters[upper*4+3];

        // reaction rates at the lower and upper pressures
        const double ln_k_lower = ln_A_lower + b_lower * ln_T - E_lower_over_R / T;
        const double ln_k_upper = ln_A_upper + b_upper * ln_T - E_upper_over_R / T;

        // interpolation
        if (p < par->plog_parameters[par->plog_seperators[j]*4+0])    // p < smallest pressure level
        {
            this->ln_k_forward[index] = ln_k_lower;
        }
        else if (par->plog_parameters[(par->plog_seperators[j+1]-1)*4+0] < p)    // p > largest pressure level
        {
            this->ln_k_forward[index] = ln_k_upper;
        }
        else
        {
            this->ln_k_forward[index] = ln_k_lower + (std::log(p) - std::log(P_lower)) / (std::log(P_upper) - std::log(P_lower)) * (ln_k_upper - ln_k_lower);
        }
    }
    
// Forward rate thresholding
    if (!cpar.enable_rate_thresholding) return;
    const double ln_correction = std::log(par->R_g * T / par->atm2Pa);  // unit correction for higher order reactions [m^3/mol]
    const double ln_uni_limit  = std::log(par->k_B * T / par->h);       // [1/s]
    const double ln_kB_T       = std::log(par->k_B) + ln_T;
    for(index_t index = 0; index < par->num_reactions; ++index)
    {
        const double ln_bi_limit_base = par->ln_bimolecular_threshold_base[index];
        const index_t m = par->reaction_order[index];
        double ln_threshold;
        if (std::isnan(ln_bi_limit_base))    // unimolecular
        {
            ln_threshold = ln_uni_limit + (m - 1) * ln_correction;
        }
        else    // bimolecular
        {
            const double T_star    = par->k_B * T / par->epsilon_AB[index];
            const double ln_T_star = ln_kB_T  - par->ln_epsilon_AB[index];
            ln_threshold = ln_bi_limit_base + 0.5 * ln_T + ln_Omega_22(T_star, ln_T_star) + (m - 2) * ln_correction;
        }
        const double ln_k = this->ln_k_forward[index];
        if (!std::isfinite(ln_k) || ln_k > ln_threshold)
        {
            this->ln_k_forward[index] = std::copysign(ln_threshold, ln_k);
        }
    }
}


void OdeFun::backward_rate(
    const double T
) 
{
    const double ln_p_over_RT = std::log(par->atm2Pa / (par->R_g * T));
    for(index_t index = 0; index < par->num_reactions; ++index)
    {
        // Equilibrium constants (K_c)
        double Delta_S = 0.0, Delta_H = 0.0;
        for (index_t k = index * par->num_max_species_per_reaction; k < (index + 1) * par->num_max_species_per_reaction; ++k)
        {
            index_t nu_index = par->nu_indexes[k];
            if (nu_index == par->invalid_index) break;

            stoich_t nu = par->nu[k];
            Delta_S += nu * this->S[nu_index];
            Delta_H += nu * this->H[nu_index];
        }

        const double ln_K_p = Delta_S / par->R_g - Delta_H / (par->R_g * T);
        double ln_K_c = ln_K_p + par->sum_nu[index] * ln_p_over_RT;
        this->ln_k_backward[index] = this->ln_k_forward[index] - ln_K_c;
    }

    // Irreversible reactions
    for(index_t j = 0; j < par->num_irreversible_reactions; ++j)
    {
        index_t index = par->irreversible_reaction_indexes[j];
        this->ln_k_backward[index] = -1e300;
    }
}


// Most stoichiometric coefficients are 0, 1, or 2
inline double pow_stoich(const double conc, const stoich_t nu)
{
    switch (nu)
    {
    case 0:
        return 1.0;
        break;
    case 1:
        return conc;
        break;
    case 2:
        return conc * conc;
        break;
    default:
        return std::pow(conc, nu);
    }
}


void OdeFun::production_rate(
    const double T,
    const double M,
    const double p,
    const double* conc
)
{
// Third body correction factors
    for (index_t j = 0; j < par->num_third_body_reactions; ++j)
    {
        double M_eff_j = 0.0;
        for (index_t k = 0; k < par->num_species; ++k)
        {
            M_eff_j += par->alpha[j*par->num_species+k] * conc[k];
        }
        this->M_eff[j] = M_eff_j;
    }
// Forward and backward rates
    this->forward_rate(T, M, p);
    this->backward_rate(T);
// Net rates
    for (index_t index = 0; index < par->num_reactions; ++index)
    {
        double forward = 1.0;
        double backward = 1.0;
        for (index_t k = index * par->num_max_species_per_reaction; k < (index + 1) * par->num_max_species_per_reaction; ++k)
        {
            index_t nu_index = par->nu_indexes[k];
            if (nu_index == par->invalid_index) break;

            forward  *= pow_stoich(conc[nu_index], par->nu_forward[k]);
            backward *= pow_stoich(conc[nu_index], par->nu_backward[k]);
        }
        const double k_forward = std::exp(this->ln_k_forward[index]);
        const double k_backward = std::exp(this->ln_k_backward[index]);
        this->net_rates[index] = k_forward * forward - k_backward * backward;
    }
// Third body reaction rates
    for (index_t j = 0; j < par->num_third_body_reactions; ++j)
    {
        if (!par->is_falloff_reaction[j])
        {
            index_t index = par->third_body_reaction_indexes[j];
            this->net_rates[index] *= this->M_eff[j];
        }
    }
// Production rates
    for (index_t k = 0; k < par->num_species; ++k)
    {
        this->omega_dot[k] = 0.0;
    }

   for (index_t index = 0; index < par->num_reactions; ++index)
   {
        for (index_t k = index * par->num_max_species_per_reaction; k < (index + 1) * par->num_max_species_per_reaction; ++k)
        {
            index_t nu_index = par->nu_indexes[k];
            if (nu_index == par->invalid_index) break;

            this->omega_dot[nu_index] += par->nu[k] * this->net_rates[index];
        }
   }
}


is_success OdeFun::operator()(
        const double t_dimless,
        const double* x_dimless,
        double* x_dimless_dot
    )
{
    if (!this->check_before_call(x_dimless))
        return false;


// Non-dimensional -> dimensional (SI)
    for (index_t k = 0; k < par->num_species+4; ++k)
    {
        this->x_dimensional[k] = x_dimless[k];
    }
    double t = t_dimless;
    cpar.dimensionalize(t, this->x_dimensional);
    double* x_dimensional_dot = x_dimless_dot;

// Common variables
    const double R = x_dimensional[0];                     // bubble radius [m]
    const double R_dot = x_dimensional[1];                 // bubble radius derivative [m/s]
    const double T = x_dimensional[2];                     // temperature [K]
    const double* conc_dimensional = x_dimensional + 3;    // molar concentrations [mol/m^3]
    double* conc_dot = x_dimensional_dot + 3;              // molar concentrations derivative [mol/m^3/s]
    double M = 0.0;                                        // sum of molar concentrations [mol/m^3]
    double C_p_avg = 0.0;                                  // average molar heat capacity at constant pressure [J/mol/K]
    double lambda_avg = 0.0;                               // average thermal conductivity [W/m/K]
    double sum_omega_dot = 0.0;                            // sum of production rates [mol/m^3/s]

// Thermodynamics
    this->thermodynamic(T);    // set C_p, H, S

// Averages
    for (index_t k = 0; k < par->num_species; ++k)
    {
        M += conc_dimensional[k];
    }
    for (index_t k = 0; k < par->num_species; ++k)
    {
        const double X_k = conc_dimensional[k] / M;
        C_p_avg += C_p[k] * X_k;
        lambda_avg += par->thermal_conductivities[k] * X_k;
    }
    const double p = this->internal_pressure(T, M, conc_dimensional);
    const double C_v_avg = C_p_avg - par->R_g;     // average molar heat capacity at constant volume [J/mol/K]

// Heat transfer
    double Q_th_dot = 0.0;
    if (cpar.enable_heat_transfer)
    {
        const double chi_avg = lambda_avg / (C_p_avg * M);
        double l_th = std::numeric_limits<double>::max();
        if (R_dot != 0.0)
        {
            l_th = std::sqrt(R * chi_avg / std::abs(R_dot));
        }
        l_th = std::min(l_th, R * std::numbers::inv_pi);
        Q_th_dot = lambda_avg * (cpar.T_inf - T) / l_th;
    }

// d/dt R
    x_dimensional_dot[0] = R_dot;

// d/dt conc
    if (cpar.enable_reactions)
    {
        this->production_rate(T, M, p, conc_dimensional);   // set omega_dot
    }
    else
    {
        std::fill(this->omega_dot, this->omega_dot + par->num_species, 0.0);
    }
    for (index_t k = 0; k < par->num_species; ++k)
    {
        conc_dot[k] = this->omega_dot[k] - conc_dimensional[k] * 3.0 * R_dot / R;
        sum_omega_dot += this->omega_dot[k];
    }

// Evaporation
    double n_net_dot = 0.0;
    double evap_energy = 0.0;
    if (cpar.enable_evaporation)
    {
        std::pair<double, double> _evap = this->evaporation(p, T, conc_dimensional[par->index_of_water]/M);
        n_net_dot = _evap.first;
        evap_energy = _evap.second;
        conc_dot[par->index_of_water] += n_net_dot * 3.0 / R;
    }

// d/dt T
    double Q_r_dot = sum_omega_dot * par->R_g * T;
    for (index_t k = 0; k < par->num_species; ++k)
    {
        Q_r_dot -= this->omega_dot[k] * this->H[k];
    }
    const double T_dot = (Q_r_dot + 3.0 / R * (-p * R_dot + Q_th_dot + evap_energy)) / (M * C_v_avg);
    const double M_dot = sum_omega_dot - 3.0 * R_dot / R * M + n_net_dot * 3.0 / R;
    const double p_dot = this->internal_pressure_derivative(T, T_dot, M, M_dot, conc_dimensional, conc_dot);
    x_dimensional_dot[2] = T_dot;

// d/dt R_dot
    const auto [P_inf, P_inf_dot] = this->excitation_pressures(t);
    const double R_dot_dot = this->bubble_dynamics(t, R, R_dot, p, p_dot, P_inf, P_inf_dot);
    x_dimensional_dot[1] = R_dot_dot;

// Dissipated energy
    if (cpar.enable_dissipated_energy)
    {
        const double V_dot = 4.0 * R * R * R_dot * std::numbers::pi;
        const double integrand_th = -(p * (1 + R_dot / cpar.c_L) + R / cpar.c_L * p_dot) * V_dot;
        const double integrand_v = 16.0 * std::numbers::pi * cpar.mu_L * (R * R_dot*R_dot + R * R * R_dot * R_dot_dot / cpar.c_L);
        const double integrand_r = 4.0 * std::numbers::pi / cpar.c_L * R * R * R_dot * (R_dot * p + p_dot * R - 0.5 * cpar.rho_L * R_dot * R_dot * R_dot - cpar.rho_L * R * R_dot * R_dot_dot);

        x_dimensional_dot[par->num_species+3] = integrand_th + integrand_v + integrand_r;
    } else {
        x_dimensional_dot[par->num_species+3] = 0.0;
    }

// Dimensional -> non-dimensionalization
    // x_dimensional_dot = x_dimless_dot
    cpar.nondimensionalize_dot(x_dimless_dot, x_dimensional);

    if (!this->check_after_call(t_dimless, x_dimless, x_dimless_dot))
        return false;
    return true;
}
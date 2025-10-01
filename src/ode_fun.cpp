#include <algorithm>
#include <sstream>
#include <iomanip>
#include <numbers>
#include <cmath>

#include "common.h"
#include "ode_fun.h"


double vapour_pressure(const double T)
{
    double T_C = T - 273.15; // [°C]
    return 611.21 * std::exp( (18.678 - T_C / 234.5) * (T_C / (257.14 + T_C)) ); // [Pa]
}


double viscosity(const double T)
{
    return 1.856e-14 * std::exp(4209.0/T + 0.04527*T - 3.376e-5*T*T); // [Pa*s]
}


OdeFun::OdeFun():
    par(nullptr),
    cpar(),
    num_species(0),
    C_p(nullptr),
    H(nullptr),
    S(nullptr),
    C_v(nullptr),
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
    if (this->C_p != nullptr)            delete[] this->C_p;
    if (this->H != nullptr)              delete[] this->H;
    if (this->S != nullptr)              delete[] this->S;
    if (this->C_v != nullptr)            delete[] this->C_v;
    if (this->M_eff != nullptr)          delete[] this->M_eff;
    if (this->ln_k_forward != nullptr)   delete[] this->ln_k_forward;
    if (this->ln_k_backward != nullptr)  delete[] this->ln_k_backward;
    if (this->net_rates != nullptr)      delete[] this->net_rates;
    if (this->omega_dot != nullptr)      delete[] this->omega_dot;

    this->C_p = nullptr;
    this->H = nullptr;
    this->S = nullptr;
    this->C_v = nullptr;
    this->M_eff = nullptr;
    this->ln_k_forward = nullptr;
    this->ln_k_backward = nullptr;
    this->net_rates = nullptr;
    this->omega_dot = nullptr;
}

is_success OdeFun::check_before_call(const double* x)
{
    // Check if iinitialization was correct
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
    if (x == nullptr)
    {
        this->cpar.error_ID = LOG_ERROR(Error::severity::error, Error::type::odefun, "State vector x is nullptr", this->cpar.ID);
        return false;
    }
    else if(!std::isfinite(x[0]) || !std::isfinite(x[2]))
    {
        std::string message = "Non finite R or T: R=" + std::to_string(x[0]) + ";   T=" + std::to_string(x[2]);
        LOG_ERROR(Error::severity::warning, Error::type::odefun, message, this->cpar.ID);  // recoverable error
        return false;
    }
    else if(x[0] < 0 || x[2] < 0)
    {
        //std::string message = "Negative R or T: R=" + std::to_string(x[0]) + ";   T=" + std::to_string(x[2]);
        //LOG_ERROR(Error::severity::warning, Error::type::odefun, message, this->cpar.ID);  // recoverable error
        return false;
    }
    
    return true;
}


is_success OdeFun::check_after_call(
    const double t,
    const double* x,
    double* dxdt
) {
    for (index_t k = 0; k < this->num_species+4; ++k)
    {
        if (!std::isfinite(dxdt[k]))
        {
            std::stringstream ss;
            if (std::isinf(dxdt[k]))
            {
                ss << "dxdt[" << k << "] is infinite";
            }
            else if (std::isnan(dxdt[k]))
            {
                ss << "dxdt[" << k << "] is NaN";
            } else
            {
                ss << "dxdt[" << k << "] is not finite";
            }
            ss << ". t = " << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10) << t;
            ss << ";    x = " << to_string((double*)x, this->num_species+4);
            ss << ";    dxdt = " << to_string((double*)dxdt, this->num_species+4);
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
    this->par = Parameters::get_parameters(cpar.mechanism);
    if (this->par == nullptr)
    {
        this->cpar.error_ID = LOG_ERROR(Error::severity::error, Error::type::odefun, "Invalid mechanism: " + std::to_string(cpar.mechanism), this->cpar.ID);
        return false;
    }
    this->num_species = this->par->num_species;
    if (this->cpar.enable_evaporation && this->par->index_of_water == this->par->invalid_index)
    {
        LOG_ERROR(Error::severity::warning, Error::type::odefun, "Water is not in the mechanism, but evaporation is enabled.", this->cpar.ID);
    }
    
    if (old_par != this->par || this->omega_dot == nullptr)
    {
        this->delete_memory();
        this->C_p           = new double[par->num_species];
        this->H             = new double[par->num_species];
        this->S             = new double[par->num_species];
        this->C_v           = new double[par->num_species];
        this->M_eff         = new double[par->num_third_bodies];
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
    const double& T_high = par->temp_range[3*par->index_of_water+1];
    const double& T_mid = par->temp_range[3*par->index_of_water+2];

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

        this->C_v_inf = par->R_erg * (a[0] + a[1]*T1 + a[2]*T2 + a[3]*T3 + a[4]*T4 - 1.0);
    }
    else
    {
    // linear extrapolation
        const double *interval_values      = &(par->interval_values[par->index_of_water*3]);       // {C_p_high, H_high, S_high}
        const double *interval_derivatives = &(par->interval_derivatives[par->index_of_water*3]);  // {dC_p_high/dT, dH_high/dT, dS_high/dT}
        this->C_v_inf = interval_values[0] + interval_derivatives[0] * (this->cpar.T_inf - T_high) - par->R_erg;
    }

    return true;
}


is_success OdeFun::initial_conditions(
    double* x
) //noexcept
{
// Equilibrium state
    const double p_E = cpar.P_amb + 2.0 * cpar.surfactant * par->sigma / cpar.R_E;   // [Pa]
    const double V_E = 4.0 / 3.0 * std::numbers::pi * cpar.R_E * cpar.R_E * cpar.R_E;    // [m^3]
    double p_gas = cpar.enable_evaporation ? p_E - cpar.P_v : p_E;   // [Pa]
    const double n_gas = p_gas * V_E / (par->R_g * cpar.T_inf);    // [mol]
    double n_H2O = cpar.enable_evaporation ? cpar.P_v * V_E / (par->R_g * cpar.T_inf) : 0.0;   // [mol]
    double c_gas = n_gas / V_E;   // [mol/m^3]
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
    x[0] = R_0;   // R_0 [m]
    x[1] = 0.0;         // R_dot_0 [m/s]
    x[2] = cpar.T_inf; // T_0 [K]
    for (index_t k = 0; k < par->num_species; ++k)
    {
        x[3+k] = 0.0;  // c_k_0 [mol/cm^3]
    }
    x[3+par->num_species] = 0.0;   // dissipated energy [J]
    if (cpar.enable_evaporation && par->index_of_water != par->invalid_index)
    {
        x[3+par->index_of_water] = 1.0e-6 * c_H2O; // c_H2O_0 [mol/cm^3]
    }
    for (index_t k = 0; k < cpar.num_initial_species; ++k)
    {
        index_t index = cpar.species[k];
        x[3+index] = 1.0e-6 * cpar.fractions[k] * c_gas;   // c_k_0 [mol/cm^3]
    }

// Errors
    if (p_gas < 0.0)
    {
        this->cpar.error_ID = LOG_ERROR(Error::severity::error, Error::type::odefun, "Negative gas pressure: " + std::to_string(p_gas), this->cpar.ID);
        return false;
    }
    return true;
}


std::pair<double, double> OdeFun::pressures(
    const double t,
    const double R,
    const double R_dot,
    const double p,
    const double p_dot
) //noexcept
{
    double p_Inf, p_Inf_dot;
    switch (cpar.excitation_type)
    {
    case Parameters::excitation::two_sinusoids:
        {
            const double& p_A1 = cpar.excitation_params[0];
            const double& p_A2 = cpar.excitation_params[1];
            const double& freq1 = cpar.excitation_params[2];
            const double& freq2 = cpar.excitation_params[3];
            const double& theta_phase = cpar.excitation_params[4];
            const double& n_cycles = cpar.excitation_params[5];

            // TODO: fix discontinuities if freq1 and freq2 are not harmonics
            if (t < 0.0 || (t > n_cycles / freq1 && t > n_cycles / freq2))
            {
                p_Inf = cpar.P_amb;
                p_Inf_dot = 0.0;
            }
            else
            {
                p_Inf = cpar.P_amb + p_A1 * sin(2.0 * std::numbers::pi * freq1 * t) + p_A2 * sin(2.0 * std::numbers::pi * freq2 * t + theta_phase);
                p_Inf_dot = 2.0 * std::numbers::pi * (p_A1 * freq1 * cos(2.0 * std::numbers::pi * freq1 * t) + p_A2 * freq2 * cos(2.0 * std::numbers::pi * freq2 * t + theta_phase));
            }
            break;
        }
    case Parameters::excitation::sin_impulse:
        {
            const double& p_A = cpar.excitation_params[0];
            const double& freq = cpar.excitation_params[1];
            const double& n_cycles = cpar.excitation_params[2];

            if (t < 0.0 || t > n_cycles / freq)
            {
                p_Inf = cpar.P_amb;
                p_Inf_dot = 0.0;
            }
            else
            {
                double insin = 2.0 * std::numbers::pi * freq;
                p_Inf = cpar.P_amb + p_A * sin(insin * t);
                p_Inf_dot = p_A * insin * cos(insin * t);
            }
            break;
        }
    case Parameters::excitation::sin_impulse_logf:
        {
            const double& p_A = cpar.excitation_params[0];
            const double& freq = pow(10.0, cpar.excitation_params[1]);
            const double& n_cycles = cpar.excitation_params[2];

            if (t < 0.0 || t > n_cycles / freq)
            {
                p_Inf = cpar.P_amb;
                p_Inf_dot = 0.0;
            }
            else
            {
                double insin = 2.0 * std::numbers::pi * freq;
                p_Inf = cpar.P_amb + p_A * sin(insin * t);
                p_Inf_dot = p_A * insin * cos(insin * t);
            }
            break;
        }
    default: // no_excitation
        {
            p_Inf = cpar.P_amb;
            p_Inf_dot = 0.0;
            break;
        }
    }
    
    const double p_L = p - (2.0 * cpar.surfactant * par->sigma + 4.0 * cpar.mu_L * R_dot) / R;
    const double p_L_dot = p_dot + (2.0 * cpar.surfactant * par->sigma * R_dot + 4.0 * cpar.mu_L * R_dot * R_dot) / (R * R); // + 4.0 * cpar.mu_L * R_ddot / R;
    const double delta = (p_L - p_Inf) / cpar.rho_L;
    const double delta_dot = (p_L_dot - p_Inf_dot) / cpar.rho_L;

    return std::make_pair(delta, delta_dot);
}


void OdeFun::thermodynamic(
    const double T
) //noexcept
{
    for (index_t k = 0; k < par->num_species; ++k)
    {
        const double& T_low = par->temp_range[3*k];
        const double& T_high = par->temp_range[3*k+1];
        const double& T_mid = par->temp_range[3*k+2];
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

        // calculate polynomials
            const double T1 = T;
            const double T2 = T1 * T1;
            const double T3 = T2 * T1;
            const double T4 = T2 * T2;
            const double T5 = T3 * T2;
        
            // Molar heat capacities at constant pressure (isobaric) [erg/mol/K]
            this->C_p[k] = par->R_erg * (a[0] + a[1]*T1 + a[2]*T2 + a[3]*T3 + a[4]*T4);
            // Enthalpies [erg/mol]
            this->H[k] = par->R_erg * (a[0]*T1 + a[1]*T2/2.0 + a[2]*T3/3.0 + a[3]*T4/4.0 + a[4]*T5/5.0 + a[5]);
            // Entropies [erg/mol/K]
            this->S[k] = par->R_erg * (a[0]*std::log(T) + a[1]*T1 + a[2]*T2/2.0 + a[3]*T3/3.0 + a[4]*T4/4.0 + a[6]);
            // Molar heat capacities at constant volume (isochoric) [erg/mol/K]
            this->C_v[k] = this->C_p[k] - par->R_erg;
        }
        else
        {
        // linear extrapolation
            const double *interval_values      = &(par->interval_values[k*3]);       // {C_p_high, H_high, S_high}
            const double *interval_derivatives = &(par->interval_derivatives[k*3]);  // {dC_p_high/dT, dH_high/dT, dS_high/dT}
            this->C_p[k] = interval_values[0] + interval_derivatives[0] * (T - T_high);
            this->H[k]   = interval_values[1] + interval_derivatives[1] * (T - T_high);
            this->S[k]   = interval_values[2] + interval_derivatives[2] * (T - T_high);
            this->C_v[k]  = this->C_p[k] - par->R_erg;
        }
    }
}


std::pair<double, double> OdeFun::evaporation(
    const double p,
    const double T,
    const double X_H2O
) //noexcept
{
// condensation and evaporation
    const double p_H2O = p * X_H2O;
    const double n_eva_dot = 1.0e3 * cpar.alfa_M * cpar.P_v / (par->W[par->index_of_water] * std::sqrt(2.0 * std::numbers::pi * par->R_v * cpar.T_inf));
    const double n_con_dot = 1.0e3 * cpar.alfa_M * p_H2O    / (par->W[par->index_of_water] * std::sqrt(2.0 * std::numbers::pi * par->R_v * T));
    const double n_net_dot = n_eva_dot - n_con_dot;
// Evaporation energy [J/mol]
    const double& C_v = this->C_v[par->index_of_water]; // Molar heat capacity of water at constant volume (isochoric) [erg/mol/K]
    double e_eva = this->C_v_inf * cpar.T_inf * 1.0e-7;
    double e_con = C_v * T * 1.0e-7;
    double evap_energy = n_eva_dot * e_eva - n_con_dot * e_con;    // [W/m^2]

    return std::make_pair(n_net_dot, evap_energy);
}


void OdeFun::forward_rate(
    const double T,
    const double M,
    const double p,
    const double ln_reaction_rate_threshold
) //noexcept
{
// Arrhenius reactions
    const double logT = std::log(T);
    for(index_t index = 0; index < par->num_reactions; ++index)
    {
        this->ln_k_forward[index] = par->ln_A[index] + par->b[index] * logT - par->E_over_R[index] / T;
    }

// Pressure dependent reactions
    index_t troe_index=0, sri_index=0;
    for(index_t j = 0; j < par->num_pressure_dependent; ++j)
    {
    // Third body reactions
        double M_eff_loc = M;
        index_t third_body_index = par->is_third_body_indexes[j];
        if (third_body_index != par->invalid_index)
        {
            M_eff_loc = this->M_eff[third_body_index];
        }
        
    // Pressure dependent formalisms
        index_t index = par->pressure_dependent_indexes[j];
        const double *reac_const = &(par->reac_const[j*3]);
        const double &ln_A0 = reac_const[0];
        const double &b_0 = reac_const[1];
        const double &E0_over_R = reac_const[2];

        const double ln_k_inf = this->ln_k_forward[index];
        const double ln_k_0 = ln_A0 + b_0 * logT - E0_over_R / T;
        const double ln_Pr = ln_k_0 + std::log(M_eff_loc) - ln_k_inf;
        const double log10_Pr = ln_Pr / std::numbers::ln10;

        double ln_F = 0.0;  // F = 1.0

        switch (par->pressure_dependent_reac_types[j])
        {
        case Parameters::reac_type::lindemann_reac:
            {
                // F = 1.0;
                break;
            }
        case Parameters::reac_type::troe_reac:
            {
                const double *troe = &(par->troe[troe_index*4]);
                const double &alfa  = troe[0];
                const double &T_xxx = troe[1];  // T***
                const double &T_x   = troe[2];  // T*
                const double &T_xx  = troe[3];  // T**

                constexpr double small_num = 1.0e-30;
                constexpr double large_num = 1.0e30;
                const double exp_negT_over_Txxx = (T_xxx <= small_num) ? 0.0 : std::exp(-T / T_xxx);
                const double exp_negT_over_Tx   = (T_x   >= large_num) ? 1.0 : std::exp(-T / T_x);
                const double exp_negTxx_over_T  = (T_xx  >= large_num) ? 0.0 : std::exp(-T_xx / T);
                const double F_cent = (1.0 - alfa) * exp_negT_over_Txxx + alfa * exp_negT_over_Tx + exp_negTxx_over_T;
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
        case Parameters::reac_type::sri_reac:   // TODO: We don't have any SRI reactions in the current mechanisms
            {
                const double *sri = &(par->sri[sri_index*5]);
                const double &a = sri[0];
                const double &b = sri[1];
                const double &c = sri[2];
                const double &d = sri[3];
                const double &e = sri[4];

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
    for(index_t j = 0; j < par->num_plog; ++j)
    {
        index_t index = par->plog_indexes[j];
        // determne indexes of the lower and upper pressures
        index_t lower = par->plog_seperators[j];
        for (index_t k = par->plog_seperators[j] + 1; k < par->plog_seperators[j+1] - 1; ++k)
        {
            const double &P_j = par->plog[k*4+0];
            if (P_j < p)
                lower = k;
            else
                break;
        }
        index_t upper = lower + 1;
        const double &P_lower        = par->plog[lower*4+0];
        const double &P_upper        = par->plog[upper*4+0];
        const double &ln_A_lower     = par->plog[lower*4+1];
        const double &ln_A_upper     = par->plog[upper*4+1];
        const double &b_lower        = par->plog[lower*4+2];
        const double &b_upper        = par->plog[upper*4+2];
        const double &E_lower_over_R = par->plog[lower*4+3];
        const double &E_upper_over_R = par->plog[upper*4+3];

        // reaction rates at the lower and upper pressures
        const double ln_k_lower = ln_A_lower + b_lower * logT - E_lower_over_R / T;
        const double ln_k_upper = ln_A_upper + b_upper * logT - E_upper_over_R / T;

        // interpolation
        if (p < par->plog[par->plog_seperators[j]*4+0])    // p < smallest pressure level
        {
            this->ln_k_forward[index] = ln_k_lower;
        }
        else if (par->plog[(par->plog_seperators[j+1]-1)*4+0] < p)    // p > largest pressure level
        {
            this->ln_k_forward[index] = ln_k_upper;
        }
        else
        {
            this->ln_k_forward[index] = ln_k_lower + (std::log(p) - std::log(P_lower)) / (std::log(P_upper) - std::log(P_lower)) * (ln_k_upper - ln_k_lower);
        }
    }

// Forward rate thresholding
    for(index_t index = 0; index < par->num_reactions; ++index)
    {
        const double ln_threshold = ln_reaction_rate_threshold + par->ln_N_A * par->reaction_order[index];
        const double ln_k_forward = this->ln_k_forward[index];
        if (!std::isfinite(ln_k_forward) || ln_k_forward > ln_threshold)
            this->ln_k_forward[index] = ln_threshold;
    }
}


void OdeFun::backward_rate(
    const double T,
    const double ln_reaction_rate_threshold
) //noexcept 
{
    for(index_t index = 0; index < par->num_reactions; ++index)
    {
        // Equilibrium constants (K_c)
        double Delta_S = 0.0, Delta_H = 0.0;
        for (index_t k = index * par->num_max_specie_per_reaction; k < (index + 1) * par->num_max_specie_per_reaction; ++k)
        {
            index_t nu_index = par->nu_indexes[k];
            if (nu_index == par->invalid_index) break;

            stoich_t nu = par->nu[k];
            Delta_S += nu * this->S[nu_index];
            Delta_H += nu * this->H[nu_index];
        }

        const double ln_K_p = Delta_S / par->R_erg - Delta_H / (par->R_erg * T);
        double ln_K_c = ln_K_p + par->sum_nu[index] * std::log(par->atm2Pa * 10.0 / (par->R_erg * T));
        double ln_k_backward = this->ln_k_forward[index] - ln_K_c;

        if (!std::isfinite(ln_K_c))
            LOG_ERROR(Error::severity::warning, Error::type::odefun, "Non finite ln_K_c in backward reaction: " + std::to_string(ln_K_c), this->cpar.ID);
        
        // Backward rate thresholding
        const double ln_threshold = ln_reaction_rate_threshold + par->ln_N_A * par->reaction_order[index] - ln_K_c;
        if (ln_k_backward > ln_threshold || !std::isfinite(ln_k_backward))
        {
            ln_k_backward = ln_threshold;
            this->ln_k_forward[index] = ln_k_backward + ln_K_c;
        }
        this->ln_k_backward[index] = ln_k_backward;
    }

    // Irreversible reactions
    for(index_t j = 0; j < par->num_irreversible; ++j)
    {
        index_t index = par->irreversible_indexes[j];
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
) //noexcept
{
// Third body correction factors
    for (index_t j = 0; j < par->num_third_bodies; ++j)
    {
        double M_eff_j = 0.0;
        for (index_t k = 0; k < par->num_species; ++k)
        {
            M_eff_j += par->alfa[j*par->num_species+k] * conc[k];
        }
        this->M_eff[j] = M_eff_j;
    }
// Forward and backward rates
    const double reaction_rate_threshold = par->k_B * T / par->h;
    const double ln_reaction_rate_threshold = std::log(reaction_rate_threshold);
    this->forward_rate(T, M, p, ln_reaction_rate_threshold);
    this->backward_rate(T, ln_reaction_rate_threshold);
// Net rates
    for (index_t index = 0; index < par->num_reactions; ++index)
    {
        double forward = 1.0;
        double backward = 1.0;
        for (index_t k = index * par->num_max_specie_per_reaction; k < (index + 1) * par->num_max_specie_per_reaction; ++k)
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
    for (index_t j = 0; j < par->num_third_bodies; ++j)
    {
        if (!par->is_pressure_dependent[j])
        {
            index_t index = par->third_body_indexes[j];
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
        for (index_t k = index * par->num_max_specie_per_reaction; k < (index + 1) * par->num_max_specie_per_reaction; ++k)
        {
            index_t nu_index = par->nu_indexes[k];
            if (nu_index == par->invalid_index) break;

            this->omega_dot[nu_index] += par->nu[k] * this->net_rates[index];
        }
   }
}


is_success OdeFun::operator()(
        const double t,
        const double* x,
        double* dxdt
    ) //noexcept
{
    if (!this->check_before_call(x))
        return false;
// Thermodynamics
    this->thermodynamic(x[2]);    // set C_p, H, S, C_v

// Common variables
    const double R = x[0];                      // bubble radius [m]
    const double R_dot = x[1];                  // bubble radius derivative [m/s]
    const double T = x[2];                      // temperature [K]
    const double* conc = x + 3;                 // molar concentrations [mol/cm^3]
    double* conc_dot = dxdt + 3;                // molar concentrations derivative [mol/cm^3/s]
    double M = 0.0;                             // sum of molar concentrations [mol/cm^3]
    double p = 0.0;                             // Partial pressure of the gases [Pa]
    double W_avg = 0.0;                         // average molecular weight [g/mol]
    double C_p_avg = 0.0;                       // average molar heat capacity at constant pressure [erg/mol/K]
    double C_v_avg = 0.0;                       // average molar heat capacity at constant volume [erg/mol/K]
    double lambda_avg = 0.0;                    // average thermal conductivity [W/m/K]
    double sum_omega_dot = 0.0;                 // sum of production rates [mol/cm^3/s]
    
    for (index_t k = 0; k < par->num_species; ++k)
    {
        M += conc[k];
    }
    p = 0.1 * M * par->R_erg * T;
    for (index_t k = 0; k < par->num_species; ++k)
    {
        const double X_k = conc[k] / M;
        W_avg += par->W[k] * X_k;
        C_p_avg += C_p[k] * X_k;
        C_v_avg += C_v[k] * X_k;
        lambda_avg += par->lambdas[k] * X_k;
    }

// Heat transfer
    double Q_th_dot = 0.0;
    if (cpar.enable_heat_transfer)
    {
        const double rho_avg = W_avg * M;
        const double chi_avg = 10.0 * lambda_avg * W_avg / (C_p_avg * rho_avg);
        double l_th = std::numeric_limits<double>::max();
        if (R_dot != 0.0)
        {
            l_th = std::sqrt(R * chi_avg / std::abs(R_dot));
        }
        l_th = std::min(l_th, R * std::numbers::inv_pi);
        Q_th_dot = lambda_avg * (cpar.T_inf - T) / l_th;
    }

// d/dt R
    dxdt[0] = R_dot;

// d/dt conc
    if (cpar.enable_reactions)
    {
        this->production_rate(T, M, p, conc);   // set omega_dot
    }
    else
    {
        std::fill(this->omega_dot, this->omega_dot + par->num_species, 0.0);
    }
    for (index_t k = 0; k < par->num_species; ++k)
    {
        conc_dot[k] = this->omega_dot[k] - conc[k] * 3.0 * R_dot / R;
        sum_omega_dot += this->omega_dot[k];
    }

// Evaporation
    double n_net_dot = 0.0;
    double evap_energy = 0.0;
    if (cpar.enable_evaporation)
    {
        std::pair<double, double> _evap = this->evaporation(p, T, conc[par->index_of_water]/M);
        n_net_dot = _evap.first;
        evap_energy = _evap.second;
        conc_dot[par->index_of_water] += 1.0e-6 * n_net_dot * 3.0 / R;
    }

// d/dt T
    double Q_r_dot = 0.0;
    for (index_t k = 0; k < par->num_species; ++k)
    {
        Q_r_dot -= this->omega_dot[k] * this->H[k];
    }
    Q_r_dot += sum_omega_dot * par->R_erg * T;
    const double T_dot = (Q_r_dot + 30.0 / R * (-p * R_dot + Q_th_dot + evap_energy)) / (M * C_v_avg);
    const double p_dot = p * (sum_omega_dot / M + T_dot / T - 3.0 * R_dot / R);
    dxdt[2] = T_dot;

// d/dt R_dot
    std::pair<double, double> _pres = this->pressures(t, R, R_dot, p, p_dot);
    const double delta     = _pres.first;
    const double delta_dot = _pres.second;

    const double nom   = (1.0 + R_dot / cpar.c_L) * delta + R / cpar.c_L * delta_dot - (1.5 - 0.5 * R_dot / cpar.c_L) * R_dot * R_dot;
    const double denom = (1.0 - R_dot / cpar.c_L) * R + 4.0 * cpar.mu_L / (cpar.c_L * cpar.rho_L);

    dxdt[1] = nom / denom;

// Dissipated energy
    if (cpar.enable_dissipated_energy)
    {
        const double V_dot = 4.0 * R * R * R_dot * std::numbers::pi;
        const double integrand_th = -(p * (1 + R_dot / cpar.c_L) + R / cpar.c_L * p_dot) * V_dot;
        const double integrand_v = 16.0 * std::numbers::pi * cpar.mu_L * (R * R_dot*R_dot + R * R * R_dot * dxdt[1] / cpar.c_L);
        const double integrand_r = 4.0 * std::numbers::pi / cpar.c_L * R * R * R_dot * (R_dot * p + p_dot * R - 0.5 * cpar.rho_L * R_dot * R_dot * R_dot - cpar.rho_L * R * R_dot * dxdt[1]);

        dxdt[par->num_species+3] = integrand_th + integrand_v + integrand_r;
    } else {
        dxdt[par->num_species+3] = 0.0;
    }

    if (!this->check_after_call(t, x, dxdt))
        return false;
    return true;
}
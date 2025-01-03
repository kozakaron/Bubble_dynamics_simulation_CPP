#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <cfenv>    // TODO: error handling

#include "common.h"
#include "ode_fun.h"


ODE::ODE():
    par(nullptr),
    cpar(nullptr),
    //x(nullptr),
    dxdt(nullptr),
    C_p(nullptr),
    H(nullptr),
    S(nullptr),
    C_v(nullptr),
    M_eff(nullptr),
    k_forward(nullptr),
    k_backward(nullptr),
    net_rates(nullptr),
    omega_dot(nullptr)
{ }


ODE::~ODE()
{
    this->delete_memory();
}


void ODE::delete_memory()
{
    if (this->cpar != nullptr)       delete this->cpar;
    //if (this->x != nullptr)          delete[] this->x;
    if (this->dxdt != nullptr)       delete[] this->dxdt;
    if (this->C_p != nullptr)        delete[] this->C_p;
    if (this->H != nullptr)          delete[] this->H;
    if (this->S != nullptr)          delete[] this->S;
    if (this->C_v != nullptr)        delete[] this->C_v;
    if (this->M_eff != nullptr)      delete[] this->M_eff;
    if (this->k_forward != nullptr)  delete[] this->k_forward;
    if (this->k_backward != nullptr) delete[] this->k_backward;
    if (this->net_rates != nullptr)  delete[] this->net_rates;
    if (this->omega_dot != nullptr)  delete[] this->omega_dot;
}


void ODE::init(const cpar_t& cpar)
{
    const Parameters *old_par = this->par;
    this->par = Parameters::get_parameters(cpar.mechanism);
    
    if (old_par != this->par || this->omega_dot == nullptr)
    {
        this->delete_memory();
        this->cpar       = new cpar_t();
        //this->x          = new double[par->num_species+4];
        this->dxdt       = new double[par->num_species+4];
        this->C_p        = new double[par->num_species];
        this->H          = new double[par->num_species];
        this->S          = new double[par->num_species];
        this->C_v        = new double[par->num_species];
        this->M_eff      = new double[par->num_third_bodies];
        this->k_forward  = new double[par->num_reactions];
        this->k_backward = new double[par->num_reactions];
        this->net_rates  = new double[par->num_reactions];
        this->omega_dot  = new double[par->num_species];
    }
    this->cpar->copy(cpar);
    
    // Evaporation calculations
    // get coefficients for T_inf
    const double *a;  // NASA coefficients (length: Parameters::NASA_order+2)
    if (this->cpar->T_inf <= par->temp_range[par->index_of_water*3+2]) // T_inf <= T_mid
    {
        a = &(par->a_low[par->index_of_water*(par->NASA_order+2)]);
    }
    else
    {
        a = &(par->a_high[par->index_of_water*(par->NASA_order+2)]);
    }
    // calculate sum
    this->C_v_inf = 0.0;
    double T_pow = 1.0;
    for (index_t n = 0; n < par->NASA_order; ++n)
    {
        this->C_v_inf += a[n] * T_pow;
        T_pow *= this->cpar->T_inf;
    }
    this->C_v_inf = par->R_erg * (this->C_v_inf - 1.0);
}


std::pair<double, double> ODE::pressures(
    const double t,
    const double R,
    const double R_dot,
    const double p,
    const double p_dot) //noexcept
{
    double p_Inf, p_Inf_dot;
    switch (cpar->excitation_type)
    {
    case Parameters::excitation::two_sinusoids:
        {
            double p_A1 = cpar->excitation_params[0];
            double p_A2 = cpar->excitation_params[1];
            double freq1 = cpar->excitation_params[2];
            double freq2 = cpar->excitation_params[3];
            double theta_phase = cpar->excitation_params[4];

            p_Inf = cpar->P_amb + p_A1 * sin(2.0 * M_PI * freq1 * t) + p_A2 * sin(2.0 * M_PI * freq2 * t + theta_phase);
            p_Inf_dot = 2.0 * M_PI * (p_A1 * freq1 * cos(2.0 * M_PI * freq1 * t) + p_A2 * freq2 * cos(2.0 * M_PI * freq2 * t + theta_phase));
            break;
        }
    case Parameters::excitation::sin_impulse:
        {
            double p_A = cpar->excitation_params[0];
            double freq = cpar->excitation_params[1];
            double n = cpar->excitation_params[2];

            if (t < 0.0 || t > n / freq)
            {
                p_Inf = cpar->P_amb;
                p_Inf_dot = 0.0;
            }
            else
            {
                double insin = 2.0 * M_PI * freq;
                p_Inf = cpar->P_amb + p_A * sin(insin * t);
                p_Inf_dot = p_A * insin * cos(insin * t);
            }
            break;
        }
    case Parameters::excitation::sin_impulse_logf:
        {
            double p_A = cpar->excitation_params[0];
            double freq = pow(10.0, cpar->excitation_params[1]);
            double n = cpar->excitation_params[2];

            if (t < 0.0 || t > n / freq)
            {
                p_Inf = cpar->P_amb;
                p_Inf_dot = 0.0;
            }
            else
            {
                double insin = 2.0 * M_PI * freq;
                p_Inf = cpar->P_amb + p_A * sin(insin * t);
                p_Inf_dot = p_A * insin * cos(insin * t);
            }
            break;
        }
    default: // no_excitation
        {
            p_Inf = cpar->P_amb;
            p_Inf_dot = 0.0;
            break;
        }
    }
    
    double p_L = p - (2.0 * cpar->surfactant * par->sigma + 4.0 * cpar->mu_L * R_dot) / R;
    double p_L_dot = p_dot + (-2.0 * cpar->surfactant * par->sigma * R_dot + 4.0 * cpar->mu_L * R_dot * R_dot) / (R * R);
    double delta = (p_L - p_Inf) / cpar->rho_L;
    double delta_dot = (p_L_dot - p_Inf_dot) / cpar->rho_L;

    return std::make_pair(delta, delta_dot);
}


void ODE::thermodynamic(
    const double T
) //noexcept
{
    const double *a;  // NASA coefficients (length: Parameters::NASA_order+2)
    for (index_t k = 0; k < par->num_species; ++k)
    {
    // get coefficients for T
        if (T <= par->temp_range[3*k+2]) // T <= T_mid
        {
            a = &(par->a_low[k*(par->NASA_order+2)]);
        }
        else
        {
            a = &(par->a_high[k*(par->NASA_order+2)]);
        }
    // calculate sums
        double C_p = 0.0, H = 0.0, S = 0.0;
        double T_pow = 1.0;
        for (index_t n = 0; n < par->NASA_order; ++n)
        {
            C_p += a[n] * T_pow;
            H += a[n] * T_pow / (n + 1);
            if (n != 0) S += a[n] * T_pow / n;
            T_pow *= T;
        }
    // calculations outside the sums
        // Molar heat capacities at constant pressure (isobaric) [erg/mol/K]
        this->C_p[k] = par->R_erg * C_p;
        // Enthalpies [erg/mol]
        this->H[k] = par->R_erg * (T * H + a[par->NASA_order]);
        // Entropies [erg/mol/K]
        this->S[k] = par->R_erg * (a[0] * std::log(T) + S + a[par->NASA_order+1]);
        // Molar heat capacities at constant volume (isochoric) [erg/mol/K]
        this->C_v[k] = this->C_p[k] - par->R_erg;
    }
}


std::pair<double, double> ODE::evaporation(
    const double p,
    const double T,
    const double X_H2O
) //noexcept
{
// condensation and evaporation
    double p_H2O = p * X_H2O;
    double n_eva_dot = 1.0e3 * cpar->alfa_M * par->P_v / (par->W[par->index_of_water] * std::sqrt(2.0 * M_PI * par->R_v * cpar->T_inf));
    double n_con_dot = 1.0e3 * cpar->alfa_M * p_H2O    / (par->W[par->index_of_water] * std::sqrt(2.0 * M_PI * par->R_v * T));
    double n_net_dot = n_eva_dot - n_con_dot;
// Molar heat capacity of water at constant volume (isochoric) [J/mol/K]
    // get coefficients for T
    const double *a;
    if (T <= par->temp_range[par->index_of_water * 3 + 2]) // T <= T_mid
    {
        a = &(par->a_low[par->index_of_water*(par->NASA_order+2)]);
    }
    else
    {
        a = &(par->a_high[par->index_of_water*(par->NASA_order+2)]);
    }
    // calculate sum
    double C_v = 0.0; double T_pow = 1.0;
    for (index_t n = 0; n < par->NASA_order; ++n)
    {
        C_v += a[n] * T_pow;
        T_pow *= T;
    }
    C_v = par->R_erg * (C_v - 1.0);
// Evaporation energy [J/mol]
    double e_eva = this->C_v_inf * cpar->T_inf * 1.0e-7;
    double e_con = C_v * T * 1.0e-7;
    double evap_energy = n_eva_dot * e_eva - n_con_dot * e_con;    // [W/m^2]

    return std::make_pair(n_net_dot, evap_energy);
}


void ODE::forward_rate(
    const double T,
    const double M,
    const double p
) //noexcept
{
// Arrhenius reactions
    for(index_t index = 0; index < par->num_reactions; ++index)
    {
        this->k_forward[index] = par->A[index] * std::pow(T, par->b[index]) * std::exp(-par->E[index] / (par->R_cal * T));
    }

// Pressure dependent reactions
    index_t troe_index=0, sri_index=0;
    for(index_t j = 0; j < par->num_pressure_dependent; ++j)
    {
        index_t index = par->pressure_dependent_indexes[j];
        double k_inf = this->k_forward[index];
        const double *reac_const = &(par->reac_const[j*3]);
        double k_0 = reac_const[0] * std::pow(T, reac_const[1]) * std::exp(-reac_const[2] / (par->R_cal * T));
        
    // Third body reactions
        double M_eff_loc = M;
        index_t third_body_index = par->is_third_body_indexes[j];
        if (third_body_index != par->invalid_index)
        {
            M_eff_loc = this->M_eff[third_body_index];
        }
        
        double P_r = k_0 / k_inf * M_eff_loc;
        double F = 1.0;

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
                double F_cent = (1.0 - troe[0]) * std::exp(-T / troe[1]) + troe[0] * std::exp(-T / troe[2]) + std::exp(-troe[3] / T);
                double logF_cent = std::log10(F_cent);
                double c2 = -0.4 - 0.67 * logF_cent;
                double n = 0.75 - 1.27 * logF_cent;
                constexpr double d = 0.14;
                double logP_r = std::log10(P_r);
                double logF = 1.0 / (1.0 + std::pow((logP_r + c2) / (n - d * (logP_r + c2)), 2)) * logF_cent;
                F = std::pow(10.0, logF);
                ++troe_index;
                break;
            }
        case Parameters::reac_type::sri_reac:   // TODO: We don't have any SRI reactions in the current mechanisms
            {
                double X = 1.0 / (1.0 + std::pow(std::log10(P_r), 2));
                const double *sri = &(par->sri[sri_index*5]);
                F = sri[3] * std::pow(sri[0] * std::exp(-sri[1] / T) + std::exp(-T / sri[2]), X) * std::pow(T, sri[4]);
                ++sri_index;
                break;
            }
        default:
            break;
        }

        this->k_forward[index] = k_inf * P_r / (1.0 + P_r) * F;
    } // pressure dependent reactions end

// PLOG reactions
    for(index_t j = 0; j < par->num_plog; ++j)
    {
        index_t index = par->plog_indexes[j];
        // determne indexes of the lower and upper pressures
        index_t lower = par->plog_seperators[j];
        for (index_t k = par->plog_seperators[j] + 1; k < par->plog_seperators[j+1] - 1; ++k)
        {
            if (par->plog[k*4+0] < p)
            {
                lower = k;
            }
        }
        index_t upper = lower + 1;

        // reaction rates at the lower and upper pressures
        double k_lower = par->plog[lower*4+1] * std::pow(T, par->plog[lower*4+2]) * std::exp(-par->plog[lower*4+3] / (par->R_cal * T));
        double k_upper = par->plog[upper*4+1] * std::pow(T, par->plog[upper*4+2]) * std::exp(-par->plog[upper*4+3] / (par->R_cal * T));

        // interpolation
        double ln_k;
        if (p < par->plog[par->plog_seperators[j]*4+0])    // p < smallest pressure level
        {
            ln_k = log(k_lower);
        }
        else if (par->plog[(par->plog_seperators[j+1]-1)*4+0] < p)    // p > largest pressure level
        {
            ln_k = log(k_upper);
        }
        else
        {
            ln_k = log(k_lower) + (log(p) - log(par->plog[lower*4+0])) / (log(par->plog[upper*4+0]) / (log(par->plog[upper*4+0]) - log(par->plog[lower*4+0]))) * (log(k_upper) - log(k_lower));
        }

        this->k_forward[index] = std::exp(ln_k);
    }
}


void ODE::backward_rate(
    const double T
) //noexcept 
{
    for(index_t index = 0; index < par->num_reactions; ++index)
    {
        double Delta_S = 0.0, Delta_H = 0.0;
        for (index_t k = index * par->num_max_specie_per_reaction; k < (index + 1) * par->num_max_specie_per_reaction; ++k)
        {
            index_t nu_index = par->nu_indexes[k];
            if (nu_index == par->invalid_index) continue;

            stoich_t nu = par->nu[k];
            Delta_S += nu * this->S[nu_index];
            Delta_H += nu * this->H[nu_index];
        }
        // TODO: fix long double overflow thing
        /*long*/ double K_p = std::exp(Delta_S / par->R_erg - Delta_H / (par->R_erg * T));
        /*long*/ double K_c = K_p * std::pow((par->atm2Pa * 10.0 / (par->R_erg * T)), par->sum_nu[index]);
        this->k_backward[index] = /*static_cast<double>*/(this->k_forward[index] / K_c);
    }
    for(index_t j = 0; j < par->num_irreversible; ++j)
    {
        index_t index = par->irreversible_indexes[j];
        this->k_backward[index] = 0.0;
    }
}


void ODE::production_rate(
    const double T,
    const double M,
    const double p,
    const double* c
) //noexcept
{
// Third body correction factors
    for (index_t j = 0; j < par->num_third_bodies; ++j)
    {
        double M_eff_j = 0.0;
        for (index_t k = 0; k < par->num_species; ++k)
        {
            M_eff_j += par->alfa[j*par->num_species+k] * c[k];
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
        for (index_t k = index * par->num_max_specie_per_reaction; k < (index + 1) * par->num_max_specie_per_reaction; ++k)
        {
            index_t nu_index = par->nu_indexes[k];
            if (nu_index == par->invalid_index) continue;

            forward  *= std::pow(c[nu_index], par->nu_forward[k]);
            backward *= std::pow(c[nu_index], par->nu_backward[k]);
        }
        this->net_rates[index] = this->k_forward[index] * forward - this->k_backward[index] * backward;
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
   std::fill(this->omega_dot, this->omega_dot + par->num_species, 0.0);
   for (index_t index = 0; index < par->num_reactions; ++index)
   {
        for (index_t k = index * par->num_max_specie_per_reaction; k < (index + 1) * par->num_max_specie_per_reaction; ++k)
        {
            index_t nu_index = par->nu_indexes[k];
            if (nu_index == par->invalid_index) continue;

            this->omega_dot[nu_index] += par->nu[k] * this->net_rates[index];
        }
   }
}


const double* ODE::operator()(
        const double t,
        const double* x
    ) //noexcept
{
// Thermodynamics
    this->thermodynamic(x[2]);    // set C_p, H, S, C_v

// Common variables
    const double R = x[0];                      // bubble radius [m]
    const double R_dot = x[1];                  // bubble radius derivative [m/s]
    const double T = x[2];                      // temperature [K]
    const double* c = x + 3;                    // molar concentrations [mol/cm^3]
    double* c_dot = this->dxdt + 3;             // molar concentrations derivative [mol/cm^3/s]
    double M = 0.0;                             // sum of molar concentrations [mol/cm^3]
    double p = 0.0;                             // Partial pressure of the gases [Pa]
    double W_avg = 0.0;                         // average molecular weight [g/mol]
    double C_p_avg = 0.0;                       // average molar heat capacity at constant pressure [erg/mol/K]
    double C_v_avg = 0.0;                       // average molar heat capacity at constant volume [erg/mol/K]
    double lambda_avg = 0.0;                    // average thermal conductivity [W/m/K]
    double sum_omega_dot = 0.0;                 // sum of production rates [mol/cm^3/s]
    
    for (index_t k = 0; k < par->num_species; ++k)
    {
        M += c[k];
    }
    p = 0.1 * M * par->R_erg * T;
    for (index_t k = 0; k < par->num_species; ++k)
    {
        const double X_k = c[k] / M;
        W_avg += par->W[k] * X_k;
        C_p_avg += C_p[k] * X_k;
        C_v_avg += C_v[k] * X_k;
        lambda_avg += par->lambdas[k] * X_k;
    }

// Heat transfer
    double Q_th_dot = 0.0;
    if (cpar->enable_heat_transfer)
    {
        const double rho_avg = W_avg * M;
        const double chi_avg = 10.0 * lambda_avg * W_avg / (C_p_avg * rho_avg);
        double l_th = (double)INFINITY;
        if (R_dot != 0.0)
        {
            l_th = std::sqrt(R * chi_avg / std::abs(R_dot));
        }
        l_th = std::min(l_th, R * M_1_PI);
        Q_th_dot = lambda_avg * (cpar->T_inf - T) / l_th;
    }

// d/dt R
    this->dxdt[0] = R_dot;

// d/dt c
    if (cpar->enable_reactions)
    {
        this->production_rate(T, M, p, c);   // set omega_dot
    }
    else
    {
        std::fill(this->omega_dot, this->omega_dot + par->num_species, 0.0);
    }
    for (index_t k = 0; k < par->num_species; ++k)
    {
        c_dot[k] = this->omega_dot[k] - c[k] * 3.0 * R_dot / R;
        sum_omega_dot += this->omega_dot[k];
    }
    std::clog << "sum_omega_dot: " << sum_omega_dot << std::endl;

// Evaporation
    double n_net_dot = 0.0;
    double evap_energy = 0.0;
    if (cpar->enable_evaporation)
    {
        std::pair<double, double> _evap = this->evaporation(p, T, c[par->index_of_water]/M);
        n_net_dot = _evap.first;
        evap_energy = _evap.second;
        c_dot[par->index_of_water] += 1.0e-6 * n_net_dot * 3.0 / R;
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
    this->dxdt[2] = T_dot;

// d/dt R_dot
    std::pair<double, double> _pres = this->pressures(t, R, R_dot, p, p_dot);
    const double delta     = _pres.first;
    const double delta_dot = _pres.second;

    const double nom   = (1.0 + R_dot / cpar->c_L) * delta + R / cpar->c_L * delta_dot - (1.5 - 0.5 * R_dot / cpar->c_L) * R_dot * R_dot;
    const double denom = (1.0 - R_dot / cpar->c_L) * R + 4.0 * cpar->mu_L / (cpar->c_L * cpar->rho_L);

    this->dxdt[1] = nom / denom;

// Dissipated energy
    if (cpar->enable_dissipated_energy)
    {
        const double V_dot = 4.0 * R * R * R_dot * M_PI;
        const double integrand_th = -(p * (1 + R_dot / cpar->c_L) + R / cpar->c_L * p_dot) * V_dot;
        const double integrand_v = 16.0 * M_PI * cpar->mu_L * (R * R_dot*R_dot + R * R * R_dot * dxdt[1] / cpar->c_L);
        const double integrand_r = 4.0 * M_PI / cpar->c_L * R * R * R_dot * (R_dot * p + p_dot * R - 0.5 * cpar->rho_L * R_dot * R_dot * R_dot - cpar->rho_L * R * R_dot * dxdt[1]);

        this->dxdt[par->num_species+3] = integrand_th + integrand_v + integrand_r;
    } else {
        this->dxdt[par->num_species+3] = 0.0;
    }

    return this->dxdt;
}

/*
TODO:
    - remove ODE::x
    - add test+benchmark for other mechanisms
    - do something with the long doubles
    - check for overflow, zero division, etc.
    - add python interface
    - OdeSolution class
    - simple RKCK45
*/
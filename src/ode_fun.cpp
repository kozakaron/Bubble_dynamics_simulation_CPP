#include "ode_fun.h"

/*________________________________control_parameters________________________________*/

ControlParameters::ControlParameters():
    ID(0),
    par(Parameters::mechanism::chemkin_ar_he),
    R_E(10.0e-06),
    ratio(1.00),
    species(nullptr),
    fractions(nullptr),
    n_species(0),
    P_amb(101325.00),
    T_inf(293.15),
    alfa_M(0.3500),
    P_v(2338.10),
    mu_L(0.001000),
    rho_L(998.20),
    c_L(1483.00),
    surfactant(1.00),
    enable_heat_transfer(true),
    enable_evaporation(true),
    enable_reactions(true),
    enable_dissipated_energy(true),
    excitation_params(nullptr)
{
    //this->target_specie = index::H2;
    //this->set_species({index::AR}, {1.0});
    this->excitation_type = Parameters::excitation::sin_impulse;
    this->set_excitation_params({-2.0e5, 30000.0, 1.0});
}

ControlParameters::~ControlParameters()
{
    if (this->species != nullptr)           delete[] this->species;
    if (this->fractions != nullptr)         delete[] this->fractions;
    if (this->excitation_params != nullptr) delete[] this->excitation_params;
}

void ControlParameters::set_species(const std::initializer_list<index_t> species_list, const std::initializer_list<double> fractions_list)
{
    // TODO: extra check: index in bound
    if (fractions_list.size() != species_list.size())
        LOG_ERROR("The number of species and fractions must be equal", this->ID);
    if (this->species != nullptr) delete[] this->species;
    if (this->fractions != nullptr) delete[] this->fractions;

    this->species = new index_t[species_list.size()];
    this->fractions = new double[fractions_list.size()];
    this->n_species = species_list.size();

    std::copy(species_list.begin(), species_list.end(), species);
    std::copy(fractions_list.begin(), fractions_list.end(), fractions);
}

void ControlParameters::set_excitation_params(const std::initializer_list<double> params_list)
{
    if (this->excitation_params != nullptr) delete[] this->excitation_params;
    if (params_list.size() != Parameters::excitation_arg_nums[this->excitation_type])
        LOG_ERROR("The number of excitation parameters must be equal to " + std::to_string(Parameters::excitation_arg_nums[this->excitation_type]), this->ID);
    this->excitation_params = new double[params_list.size()];
    std::copy(params_list.begin(), params_list.end(), excitation_params);
}

void ControlParameters::copy(const ControlParameters& cpar)
{
    this->ID = cpar.ID;
    this->par = cpar.par;
    this->R_E = cpar.R_E;
    this->ratio = cpar.ratio;
    this->n_species = cpar.n_species;
    this->P_amb = cpar.P_amb;
    this->T_inf = cpar.T_inf;
    this->alfa_M = cpar.alfa_M;
    this->P_v = cpar.P_v;
    this->mu_L = cpar.mu_L;
    this->rho_L = cpar.rho_L;
    this->c_L = cpar.c_L;
    this->surfactant = cpar.surfactant;
    this->enable_heat_transfer = cpar.enable_heat_transfer;
    this->enable_evaporation = cpar.enable_evaporation;
    this->enable_reactions = cpar.enable_reactions;
    this->enable_dissipated_energy = cpar.enable_dissipated_energy;
    this->target_specie = cpar.target_specie;

    if (this->species != nullptr) delete[] this->species;
    if (this->fractions != nullptr) delete[] this->fractions;
    if (this->excitation_params != nullptr) delete[] this->excitation_params;

    this->species = new index_t[cpar.n_species];
    this->fractions = new double[cpar.n_species];
    this->excitation_params = new double[Parameters::excitation_arg_nums[cpar.excitation_type]];

    std::copy(cpar.species, cpar.species + cpar.n_species, this->species);
    std::copy(cpar.fractions, cpar.fractions + cpar.n_species, this->fractions);
    std::copy(cpar.excitation_params, cpar.excitation_params + Parameters::excitation_arg_nums[cpar.excitation_type], this->excitation_params);
    this->excitation_type = cpar.excitation_type;
}



/*________________________________ode_function________________________________*/

ODE::ODE():
    par(nullptr),
    cpar(nullptr),
    x(nullptr),
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
    if (this->x != nullptr)          delete[] this->x;
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
    this->delete_memory();
    this->cpar       = new cpar_t();
    this->par = Parameters::get_parameters(cpar.par);
    this->cpar->copy(cpar);
    
    this->x          = new double[par->num_species+4];
    this->C_p        = new double[par->num_species];
    this->H          = new double[par->num_species];
    this->S          = new double[par->num_species];
    this->C_v        = new double[par->num_species];
    this->M_eff      = new double[par->num_third_bodies];
    this->k_forward  = new double[par->num_reactions];
    this->k_backward = new double[par->num_reactions];
    this->net_rates  = new double[par->num_reactions];
    this->omega_dot  = new double[par->num_species];
    // TODO: evaporation calculations
}


std::pair<double, double> ODE::pressures(
    const double t,
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
        }
    default: // no_excitation
        {
            p_Inf = cpar->P_amb;
            p_Inf_dot = 0.0;
            break;
        }
    }
    
    double p_L = p - (2.0 * cpar->surfactant * par->sigma + 4.0 * cpar->mu_L * this->x[1]) / this->x[0];
    double p_L_dot = p_dot + (-2.0 * cpar->surfactant * par->sigma * this->x[1] + 4.0 * cpar->mu_L * this->x[1] * this->x[1]) / (this->x[0] * this->x[0]);
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
    // TODO: remember C_V_inf
    // get coefficients for T_inf
    if (cpar->T_inf <= par->temp_range[par->index_of_water*3+2]) // T_inf <= T_mid
    {
        a = &(par->a_low[par->index_of_water*(par->NASA_order+2)]);
    }
    else
    {
        a = &(par->a_high[par->index_of_water*(par->NASA_order+2)]);
    }
    // calculate sum
    double C_v_inf = 0.0; T_pow = 1.0;
    for (index_t n = 0; n < par->NASA_order; ++n)
    {
        C_v_inf += a[n] * T_pow;
        T_pow *= cpar->T_inf;
    }
    C_v_inf = par->R_erg * (C_v_inf - 1.0);
// Evaporation energy [J/mol]
    double e_eva = C_v_inf * cpar->T_inf * 1.0e-7;
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
        double k_0 = par->reac_const[j*3+0] * std::pow(T, par->reac_const[j*3+1]) * std::exp(-par->reac_const[j*3+2] / (par->R_cal * T));
        
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
                double F_cent = (1.0 - par->troe[troe_index*4+0]) * std::exp(-T / par->troe[troe_index*4+1]) + par->troe[troe_index*4+0] * std::exp(-T / par->troe[troe_index*4+2]) + std::exp(-par->troe[troe_index*4+3] / T);
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
                F = par->sri[sri_index*5+3] * std::pow(par->sri[sri_index*5+0] * std::exp(-par->sri[sri_index*5+1] / T) + std::exp(-T / par->sri[sri_index*5+2]), X) * std::pow(T, par->sri[sri_index*5+4]);
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
        char sum_nu_i = 0.0;
        for(index_t k = 0; k < par->num_max_specie_per_reaction; ++k)
        {
            index_t nu_index = par->nu_indexes[index*par->num_max_specie_per_reaction+k];
            if (nu_index == par->invalid_index) continue;

            stoich_t nu = par->nu[(index*3+2)*par->num_max_specie_per_reaction+k];
            Delta_S += nu * this->S[nu_index];
            Delta_H += nu * this->H[nu_index];
            sum_nu_i += nu;
        }
        // TODO: fix long double overflow thing
        /*long*/ double K_p = std::exp(Delta_S / par->R_erg - Delta_H / (par->R_erg * T));
        /*long*/ double K_c = K_p * std::pow((par->atm2Pa * 10.0 / (par->R_erg * T)), sum_nu_i);
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
    const double p
) //noexcept
{
// Third body correction factors
    double *c = this->x + 3;
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
        for (index_t k = 0; k < par->num_max_specie_per_reaction; ++k)
        {
            index_t nu_index = par->nu_indexes[index*par->num_max_specie_per_reaction+k];
            if (nu_index == par->invalid_index) continue;

            index_t idx_forward  = (index*3+0)*par->num_max_specie_per_reaction+k;
            index_t idx_backward = (index*3+1)*par->num_max_specie_per_reaction+k;
            forward  *= std::pow(c[nu_index], par->nu[idx_forward]);
            backward *= std::pow(c[nu_index], par->nu[idx_backward]);
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
        for(index_t k = 0; k < par->num_max_specie_per_reaction; ++k)
        {
            index_t nu_index = par->nu_indexes[index*par->num_max_specie_per_reaction+k];
            if (nu_index == par->invalid_index) continue;

            index_t idx_nu = (index*3+2)*par->num_max_specie_per_reaction+k;
            this->omega_dot[nu_index] += par->nu[idx_nu] * this->net_rates[index];
        }
   }
}

/*
TODO:
    - implement and test check_cpar()
    - do something with the long doubles
*/
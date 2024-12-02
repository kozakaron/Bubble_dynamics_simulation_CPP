#include "ode_fun.h"

/*________________________________control_parameters________________________________*/

ControlParameters::ControlParameters():
    ID(0),
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
#if defined CHEMKIN_AR_HE_H
    this->target_specie = par::index::H2;
    this->set_species({par::index::AR}, {1.0});
#elif defined CHEMKIN_OTOMO2018_H || defined CHEMKIN_OTOMO2018_WITHOUT_O_H
    this->target_specie = par::index::NH3;
    this->set_species({par::index::H2, par::index::N2}, {0.75, 0.25});
#else
    #error Unknown mechanism. See ControlParameters::ControlParameters() in ode_fun.cpp.
#endif
    this->excitation_type = par::excitation::sin_impulse;
    this->set_excitation_params({-2.0e5, 30000.0, 1.0});
}

ControlParameters::~ControlParameters()
{
    if (this->species != nullptr)           delete[] this->species;
    if (this->fractions != nullptr)         delete[] this->fractions;
    if (this->excitation_params != nullptr) delete[] this->excitation_params;
}

void ControlParameters::set_species(const std::initializer_list<par::index> species_list, const std::initializer_list<double> fractions_list)
{
    if (fractions_list.size() != species_list.size())
        LOG_ERROR("The number of species and fractions must be equal", this->ID);
    if (this->species != nullptr) delete[] this->species;
    if (this->fractions != nullptr) delete[] this->fractions;

    this->species = new par::index[species_list.size()];
    this->fractions = new double[fractions_list.size()];
    this->n_species = species_list.size();

    std::copy(species_list.begin(), species_list.end(), species);
    std::copy(fractions_list.begin(), fractions_list.end(), fractions);
}

void ControlParameters::set_excitation_params(const std::initializer_list<double> params_list)
{
    if (this->excitation_params != nullptr) delete[] this->excitation_params;
    if (params_list.size() != par::excitation_arg_nums[this->excitation_type])
        LOG_ERROR("The number of excitation parameters must be equal to " + std::to_string(par::excitation_arg_nums[this->excitation_type]), this->ID);
    this->excitation_params = new double[params_list.size()];
    std::copy(params_list.begin(), params_list.end(), excitation_params);
}

void ControlParameters::copy(const ControlParameters& cpar)
{
    this->ID = cpar.ID;
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

    this->species = new par::index[cpar.n_species];
    this->fractions = new double[cpar.n_species];
    this->excitation_params = new double[par::excitation_arg_nums[cpar.excitation_type]];

    std::copy(cpar.species, cpar.species + cpar.n_species, this->species);
    std::copy(cpar.fractions, cpar.fractions + cpar.n_species, this->fractions);
    std::copy(cpar.excitation_params, cpar.excitation_params + par::excitation_arg_nums[cpar.excitation_type], this->excitation_params);
    this->excitation_type = cpar.excitation_type;
}



/*________________________________ode_function________________________________*/

ODE::ODE()
{
    this->cpar       = new cpar_t();
    this->x          = new double[par::num_species+4];
    this->C_p        = new double[par::num_species];
    this->H          = new double[par::num_species];
    this->S          = new double[par::num_species];
    this->C_v        = new double[par::num_species];
    this->M_eff      = new double[par::num_third_bodies];
    this->k_forward  = new double[par::num_reactions];
    this->k_backward = new double[par::num_reactions];
}
ODE::~ODE()
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
}


std::pair<double, double> ODE::pressures(
    const double t,
    const double p,
    const double p_dot) //noexcept
{
    double p_Inf, p_Inf_dot;
    switch (this->cpar->excitation_type)
    {
    case par::excitation::two_sinusoids:
        {
            // TODO
            p_Inf = this->cpar->P_amb;
            p_Inf_dot = 0.0;
            break;
        }
    case par::excitation::sin_impulse:
        {
            double p_A = this->cpar->excitation_params[0];
            double freq = this->cpar->excitation_params[1];
            double n = this->cpar->excitation_params[2];
            if (t < 0.0 || t > n / freq)
            {
                p_Inf = this->cpar->P_amb;
                p_Inf_dot = 0.0;
            }
            else
            {
                double insin = 2.0 * M_PI * freq;
                p_Inf = this->cpar->P_amb + p_A * sin(insin * t);
                p_Inf_dot = p_A * insin * cos(insin * t);
            }
            break;
        }
    default: // no_excitation
        {
            p_Inf = this->cpar->P_amb;
            p_Inf_dot = 0.0;
            break;
        }
    }
    
    double p_L = p - (2.0 * this->cpar->surfactant * par::sigma + 4.0 * this->cpar->mu_L * this->x[1]) / this->x[0];
    double p_L_dot = p_dot + (-2.0 * this->cpar->surfactant * par::sigma * this->x[1] + 4.0 * this->cpar->mu_L * this->x[1] * this->x[1]) / (this->x[0] * this->x[0]);
    double delta = (p_L - p_Inf) / this->cpar->rho_L;
    double delta_dot = (p_L_dot - p_Inf_dot) / this->cpar->rho_L;

    return std::make_pair(delta, delta_dot);
}


void ODE::thermodynamic(
    const double T
) //noexcept
{
    const double *a;  // NASA coefficients (length: par::NASA_order+2)
    for (size_t k = 0; k < par::num_species; ++k)
    {
    // get coefficients for T
        if (T <= par::temp_range[k][2]) // T <= T_mid
        {
            a = par::a_low[k].data();
        }
        else
        {
            a = par::a_high[k].data();
        }
    // calculate sums
        double C_p = 0.0, H = 0.0, S = 0.0;
        double T_pow = 1.0;
        for (size_t n = 0; n < par::NASA_order; ++n)
        {
            C_p += a[n] * T_pow;
            H += a[n] * T_pow / (n + 1);
            S += (n != 0) * a[n] * T_pow / n;
            T_pow *= T;
        }
    // calculations outside the sums
        // Molar heat capacities at constant pressure (isobaric) [erg/mol/K]
        this->C_p[k] = par::R_erg * C_p;
        // Enthalpies [erg/mol]
        this->H[k] = par::R_erg * (T * H + a[par::NASA_order]);
        // Entropies [erg/mol/K]
        this->S[k] = par::R_erg * (a[0] * log(T) + S + a[par::NASA_order+1]);
        // Molar heat capacities at constant volume (isochoric) [erg/mol/K]
        this->C_v[k] = this->C_p[k] - par::R_erg;
    }
}


std::pair<double, double> ODE::evaporation(
    const double p,
    const double T
) //noexcept
{
    // TODO
    (void)p; (void)T;
    return std::make_pair<double, double>(0.0, 0.0);
}


void ODE::forward_rate(
    const double T,
    const double M,
    const double p
) //noexcept
{
// Arrhenius reactions
    for(size_t index = 0; index < par::num_reactions; ++index)
    {
        this->k_forward[index] = par::A[index] * std::pow(T, par::b[index]) * std::exp(-par::E[index] / (par::R_cal * T));
    }

// Pressure dependent reactions
    size_t troe_index=0, sri_index=0;
    for(size_t j = 0; j < par::num_pressure_dependent; ++j)
    {
        size_t index = par::pressure_dependent_indexes[j];
        double k_inf = this->k_forward[index];
        double k_0 = par::reac_const[j][0] * std::pow(T, par::reac_const[j][1]) * std::exp(-par::reac_const[j][2] / (par::R_cal * T));
        
    // Third body reactions
        double M_eff_loc = M;
        for (size_t i = 0; i < par::num_third_bodies; ++i)
        {
            if (par::third_body_indexes[i] == index)
            {
                M_eff_loc = this->M_eff[i];
                break;
            }
        }
        
        double P_r = k_0 / k_inf * M_eff_loc;
        double F = 1.0;

        switch (par::pressure_dependent_reac_types[j])
        {
        case par::reac_type::LINDEMANN:
            {
                // F = 1.0;
                break;
            }
        case par::reac_type::TROE:
            {
                double F_cent = (1.0 - par::troe[troe_index][0]) * std::exp(-T / par::troe[troe_index][1]) + par::troe[troe_index][0] * std::exp(-T / par::troe[troe_index][2]) + std::exp(-par::troe[troe_index][3] / T);
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
        case par::reac_type::SRI:   // TODO: We don't have any SRI reactions in the current mechanisms
            {
                double X = 1.0 / (1.0 + std::pow(std::log10(P_r), 2));
                F = par::sri[sri_index][3] * std::pow(par::sri[sri_index][0] * std::exp(-par::sri[sri_index][1] / T) + std::exp(-T / par::sri[sri_index][2]), X) * std::pow(T, par::sri[sri_index][4]);
                ++sri_index;
                break;
            }
        default:
            break;
        }

        this->k_forward[index] = k_inf * P_r / (1.0 + P_r) * F;
    } // pressure dependent reactions end

// PLOG reactions
    for(size_t j = 0; j < par::num_plog; ++j)
    {
        size_t index = par::plog_indexes[j];
        // determne indexes of the lower and upper pressures
        size_t lower = par::plog_seperators[j];
        for (size_t k = par::plog_seperators[j] + 1; k < par::plog_seperators[j+1] - 1; ++k)
        {
            if (par::plog[k][0] < p)
            {
                lower = k;
            }
        }
        size_t upper = lower + 1;

        // reaction rates at the lower and upper pressures
        double k_lower = par::plog[lower][1] * std::pow(T, par::plog[lower][2]) * std::exp(-par::plog[lower][3] / (par::R_cal * T));
        double k_upper = par::plog[upper][1] * std::pow(T, par::plog[upper][2]) * std::exp(-par::plog[upper][3] / (par::R_cal * T));

        // interpolation
        double ln_k;
        if (p < par::plog[par::plog_seperators[j]][0])    // p < smallest pressure level
        {
            ln_k = log(k_lower);
        }
        else if (par::plog[par::plog_seperators[j+1] - 1][0] < p)    // p > largest pressure level
        {
            ln_k = log(k_upper);
        }
        else
        {
            ln_k = log(k_lower) + (log(p) - log(par::plog[lower][0])) / (log(par::plog[upper][0]) / (log(par::plog[upper][0]) - log(par::plog[lower][0]))) * (log(k_upper) - log(k_lower));
        }

        this->k_forward[index] = std::exp(ln_k);
    }
}


void ODE::backward_rate(
    const double T
) //noexcept 
{
    for(size_t index = 0; index < par::num_reactions; ++index)
    {
        double Delta_S = 0.0, Delta_H = 0.0;
        size_t sum_nu_i = 0.0;
        for(size_t k = 0; k < par::num_species; ++k)
        {
            Delta_S += par::nu[index][k] * this->S[k];
            Delta_H += par::nu[index][k] * this->H[k];
            sum_nu_i += par::nu[index][k];
        }
        long double K_p = std::exp(Delta_S / par::R_erg - Delta_H / (par::R_erg * T));
        long double K_c = K_p * std::pow((par::atm2Pa * 10.0 / (par::R_erg * T)), sum_nu_i);
        this->k_backward[index] = static_cast<double>(this->k_forward[index] / K_c);
    }
    for(size_t j = 0; j < par::num_irreversible; ++j)
    {
        size_t index = par::irreversible_indexes[j];
        this->k_backward[index] = 0.0;
    }
}
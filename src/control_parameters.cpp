#include <sstream>
#include <numeric>

#include "control_parameters.h"


ControlParameters::ControlParameters():
    ID(0),
    mechanism(Parameters::mechanism::chemkin_otomo2018),
    error_ID(ErrorHandler::no_error),
    R_E(10.0e-06),
    species(nullptr),
    fractions(nullptr),
    num_initial_species(0),
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
    const Parameters* par = Parameters::get_parameters(Parameters::mechanism::chemkin_otomo2018);
    this->target_specie = par->get_species("NH3");
    this->set_species({par->get_species("H2"), par->get_species("N2")}, {0.75, 0.25});
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
    const Parameters* par = Parameters::get_parameters(this->mechanism);
    if (par == nullptr)
    {
        this->error_ID = LOG_ERROR("Invalid mechanism: " + std::to_string(this->mechanism), this->ID);
        return;
    }
    for (const auto& species: species_list)
        if (species == par->invalid_index)
        {
            this->error_ID = LOG_ERROR("Invalid species (" + std::to_string(species) + ") for mechanism " + par->model, this->ID);
            return;
        }
        else if (species >= par->num_species)
        {
            this->error_ID = LOG_ERROR("Species index " + std::to_string(species) + " out of bound for mechanism " + par->model, this->ID);
            return;
        }
    if (fractions_list.size() != species_list.size())
    {
        std::stringstream ss;
        ss << "The number of species and fractions must be equal: species_list=[ ";
        for (const auto& species: species_list) ss << species << " ";
        ss << "], fractions_list=[ ";
        for (const auto& fraction: fractions_list) ss << fraction << " ";
        ss << "]";
        this->error_ID = LOG_ERROR(ss.str(), this->ID);
        return;
    }
    double sum_fraction = std::accumulate(fractions_list.begin(), fractions_list.end(), 0.0);
    if (std::abs(sum_fraction - 1.0) > 1.0e-10)
    {
        std::stringstream ss;
        ss << "The sum of fractions must be equal to 1.0, instead it is " << sum_fraction << ": fractions_list=[ ";
        for (const auto& fraction: fractions_list) ss << fraction << " ";
        ss << "]";
        this->error_ID = LOG_ERROR(ss.str(), this->ID);
        return;
    }

    if (this->species != nullptr) delete[] this->species;
    if (this->fractions != nullptr) delete[] this->fractions;
    this->species = new index_t[species_list.size()];
    this->fractions = new double[fractions_list.size()];
    this->num_initial_species = species_list.size();

    std::copy(species_list.begin(), species_list.end(), species);
    std::copy(fractions_list.begin(), fractions_list.end(), fractions);
}


void ControlParameters::set_excitation_params(const std::initializer_list<double> params_list)
{
    if (params_list.size() != Parameters::excitation_arg_nums[this->excitation_type])
    {
        std::stringstream ss;
        ss << "The number of excitation parameters must be equal to " << Parameters::excitation_arg_nums[this->excitation_type] << ": params_list=[ ";
        for (const auto& param: params_list) ss << param << " ";
        ss << "]";
        this->error_ID = LOG_ERROR(ss.str(), this->ID);
        return;
    }
    if (this->excitation_params != nullptr) delete[] this->excitation_params;
    this->excitation_params = new double[params_list.size()];
    std::copy(params_list.begin(), params_list.end(), excitation_params);
}


void ControlParameters::copy(const ControlParameters& cpar)
{
    this->ID = cpar.ID;
    this->mechanism = cpar.mechanism;
    this->R_E = cpar.R_E;
    this->num_initial_species = cpar.num_initial_species;
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

    this->species = new index_t[cpar.num_initial_species];
    this->fractions = new double[cpar.num_initial_species];
    this->excitation_params = new double[Parameters::excitation_arg_nums[cpar.excitation_type]];

    std::copy(cpar.species, cpar.species + cpar.num_initial_species, this->species);
    std::copy(cpar.fractions, cpar.fractions + cpar.num_initial_species, this->fractions);
    std::copy(cpar.excitation_params, cpar.excitation_params + Parameters::excitation_arg_nums[cpar.excitation_type], this->excitation_params);
    this->excitation_type = cpar.excitation_type;
}
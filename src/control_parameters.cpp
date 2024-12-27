#include "control_parameters.h"


ControlParameters::ControlParameters():
    ID(0),
    mechanism(Parameters::mechanism::chemkin_otomo2018),
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
    for (const auto& species: species_list)
        if (species >= par->num_species)
            LOG_ERROR("Species index " + std::to_string(species) + " out of bound", this->ID);
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
    this->mechanism = cpar.mechanism;
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
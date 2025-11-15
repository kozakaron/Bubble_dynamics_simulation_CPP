#include <string>
#include <sstream>
#include <fstream>
#include <numeric>
#include <numbers>
#include <iomanip>
#include <stdint.h>

#include "nlohmann/json.hpp"
#include "common.h"
#include "parameters.h"
#include "control_parameters.h"

using ordered_json = nlohmann::ordered_json;


ControlParameters::ControlParameters()
{
    this->par = nullptr;
    this->excitation_type = Parameters::excitation::no_excitation;
}


ControlParameters::ControlParameters(const Builder& builder)
{
    this->init(builder);
}


// Helper function to get value from JSON with type checking
template<typename T>
T get_value(const ordered_json& j, const std::string& key, const T& default_value)
{
    constexpr bool is_floating_point = std::is_same_v<T, double> || std::is_same_v<T, float>;
    constexpr bool is_integer = std::is_same_v<T, int8_t> || std::is_same_v<T, int16_t> || std::is_same_v<T, int32_t> || std::is_same_v<T, int64_t> ||
                                std::is_same_v<T, uint8_t> || std::is_same_v<T, uint16_t> || std::is_same_v<T, uint32_t> || std::is_same_v<T, uint64_t>;
    constexpr bool is_string = std::is_same_v<T, std::string>;
    constexpr bool is_bool = std::is_same_v<T, bool>;
    constexpr bool is_float_vector = std::is_same_v<T, std::vector<double>> || std::is_same_v<T, std::vector<float>>;
    constexpr bool is_string_vector = std::is_same_v<T, std::vector<std::string>>;

    // check existence
    if (!j.contains(key))
    {
        LOG_ERROR(
            Error::severity::warning, Error::type::preprocess,
            "Key \"" + key + "\" not found in JSON object. Using default value. ");
        return default_value;
    }

    // check type
    std::string message = "";
    if (is_floating_point && !j.at(key).is_number_float())
    {
        message = "Expected floating point number for key \"" + key + "\", instead found " + j.at(key).dump() + ". Using default value.";
    }
    else if (is_integer && !j.at(key).is_number_integer())
    {
        message = "Expected integer number for key \"" + key + "\", instead found " + j.at(key).dump() + ". Using default value.";
    }
    else if (is_string && !j.at(key).is_string())
    {
        message = "Expected string for key \"" + key + "\", instead found " + j.at(key).dump() + ". Using default value.";
    }
    else if (is_bool && !j.at(key).is_boolean())
    {
        message = "Expected boolean for key \"" + key + "\", instead found " + j.at(key).dump() + ". Using default value.";
    }
    else if ((is_float_vector || is_string_vector) && !j.at(key).is_array())
    {
        message = "Expected JSON array for key \"" + key + "\", instead found " + j.at(key).dump() + ". Using default value.";
    }

    if (!message.empty())
    {
        LOG_ERROR(Error::severity::warning, Error::type::preprocess,message);
        return default_value;
    }

    // check element types for vectors
    if (is_float_vector || is_string_vector)
    {
        for (const auto& element : j.at(key))
        {
            if (is_float_vector && !element.is_number_float())
            {
                LOG_ERROR(
                    Error::severity::warning, Error::type::preprocess,
                    "Expected floating point number in array for key \"" + key + "\", instead found " + element.dump() + ". Using default value."
                );
                return default_value;
            }
            else if (is_string_vector && !element.is_string())
            {
                LOG_ERROR(
                    Error::severity::warning, Error::type::preprocess,
                    "Expected string in array for key \"" + key + "\", instead found " + element.dump() + ". Using default value."
                );
                return default_value;
            }
        }
    }

    // get element
    return j.at(key).get<T>();
}


ControlParameters::ControlParameters(const ordered_json& j)
{
    auto builder = ControlParameters::Builder{};
    try
    {
        builder.ID =                        get_value<size_t>                   (j, "ID",                       builder.ID);
        builder.mechanism =                 get_value<std::string>              (j, "mechanism",                builder.mechanism);
        builder.R_E =                       get_value<double>                   (j, "R_E",                      builder.R_E);
        builder.ratio =                     get_value<double>                   (j, "ratio",                    builder.ratio);
        builder.species =                   get_value<std::vector<std::string>> (j, "species",                  builder.species);
        builder.fractions =                 get_value<std::vector<double>>      (j, "fractions",                builder.fractions);
        builder.P_amb =                     get_value<double>                   (j, "P_amb",                    builder.P_amb);
        builder.T_inf =                     get_value<double>                   (j, "T_inf",                    builder.T_inf);
        builder.alpha_M =                   get_value<double>                   (j, "alpha_M",                  builder.alpha_M);
        builder.P_v =                       get_value<double>                   (j, "P_v",                      builder.P_v);
        builder.mu_L =                      get_value<double>                   (j, "mu_L",                     builder.mu_L);
        builder.rho_L =                     get_value<double>                   (j, "rho_L",                    builder.rho_L);
        builder.c_L =                       get_value<double>                   (j, "c_L",                      builder.c_L);
        builder.surfactant =                get_value<double>                   (j, "surfactant",               builder.surfactant);
        builder.enable_heat_transfer =      get_value<bool>                     (j, "enable_heat_transfer",     builder.enable_heat_transfer);
        builder.enable_evaporation =        get_value<bool>                     (j, "enable_evaporation",       builder.enable_evaporation);
        builder.enable_reactions =          get_value<bool>                     (j, "enable_reactions",         builder.enable_reactions);
        builder.enable_dissipated_energy =  get_value<bool>                     (j, "enable_dissipated_energy", builder.enable_dissipated_energy);
        builder.enable_van_der_waals =      get_value<bool>                     (j, "enable_van_der_waals",     builder.enable_van_der_waals);
        builder.enable_rate_thresholding =  get_value<bool>                     (j, "enable_rate_thresholding", builder.enable_rate_thresholding);
        builder.target_specie =             get_value<std::string>              (j, "target_specie",            builder.target_specie);
        builder.excitation_type = Parameters::string_to_excitation(
                                            get_value<std::string>              (j, "excitation_type",          Parameters::excitation_names.at(builder.excitation_type))
        );
        builder.excitation_params =         get_value<std::vector<double>>      (j, "excitation_params",        builder.excitation_params);
        builder.excitation_cycles =         get_value<double>                   (j, "excitation_cycles",        builder.excitation_cycles);
        builder.ramp_up_cycles =            get_value<double>                   (j, "ramp_up_cycles",           builder.ramp_up_cycles);
    }
    catch(const std::exception& e)
    {
        builder.error_ID = LOG_ERROR(
            Error::severity::error,
            Error::type::preprocess,
            "Error parsing JSON file: " + std::string(e.what())
        );
    }

    this->init(builder);
}


ControlParameters::ControlParameters(const std::string& json_path)
{
    // Open JSON file
    std::ifstream input_file(json_path);
    if (!input_file.is_open())
    {
        *this = ControlParameters();
        this->error_ID = LOG_ERROR(
            Error::severity::error,
            Error::type::preprocess,
            "Could not open JSON file: " + json_path
        );
        return;
    }

    // Read JSON data
    ordered_json j;
    try
    {
        j = ordered_json::parse(input_file);
    }
    catch (const std::exception& e)
    {
        *this = ControlParameters();
        this->error_ID = LOG_ERROR(
            Error::severity::error,
            Error::type::preprocess,
            "Error parsing JSON file: " + std::string(e.what())
        );
        return;
    }
    input_file.close();

    // Parse JSON data
    if (!j.contains("cpar"))
    {
        *this = ControlParameters();
        this->error_ID = LOG_ERROR(
            Error::severity::error,
            Error::type::preprocess,
            "JSON file does not contain 'cpar' key: " + j.dump()
        );
        return;
    }

    *this = ControlParameters(j.at("cpar"));
}


void ControlParameters::init(const ControlParameters::Builder& builder)
{
    this->ID = builder.ID;
    this->error_ID = builder.error_ID;
    this->set_mechanism(builder.mechanism);
    this->R_E = builder.R_E;
    this->ratio = builder.ratio;
    this->P_amb = builder.P_amb;
    this->T_inf = builder.T_inf;
    this->alpha_M = builder.alpha_M;
    this->P_v = builder.P_v;
    this->mu_L = builder.mu_L;
    this->rho_L = builder.rho_L;
    this->c_L = builder.c_L;
    this->surfactant = builder.surfactant;
    this->enable_heat_transfer = builder.enable_heat_transfer;
    this->enable_evaporation = builder.enable_evaporation;
    this->enable_reactions = builder.enable_reactions;
    this->enable_dissipated_energy = builder.enable_dissipated_energy;
    this->enable_van_der_waals = builder.enable_van_der_waals;
    this->enable_rate_thresholding = builder.enable_rate_thresholding;
    this->target_specie = par->get_species(builder.target_specie);
    this->excitation_type = builder.excitation_type;
    this->excitation_cycles = builder.excitation_cycles;
    this->ramp_up_cycles = builder.ramp_up_cycles;

    if (this->target_specie == par->invalid_index)
    {
        std::string message = "Invalid target specie (" + builder.target_specie + ") for mechanism " + par->mechanism_name;
        this->error_ID = LOG_ERROR(Error::severity::error, Error::type::preprocess, message, this->ID);
        return;
    }
    this->num_initial_species = 0;
    this->set_species(std::vector<std::string>(builder.species), std::vector<double>(builder.fractions));
    this->set_excitation_params(std::vector<double>(builder.excitation_params));

    // Set reference values
    const double p_E = this->P_amb + 2.0 * this->surfactant * par->sigma / this->R_E;   // [Pa]
    const double V_E = 4.0 / 3.0 * std::numbers::pi * this->R_E * this->R_E * this->R_E;    // [m^3]
    const double n_gas = p_E * V_E / (par->R_g * this->T_inf);    // [mol]
    double M = n_gas / V_E;   // [mol/m^3]
    (void)M;

    this->R_ref = this->R_E;
    this->T_ref = this->T_inf;
    // J = kg*m^2/s^2
    this->E_diss_ref = 1.0;
    this->c_ref = n_gas;  // TODO: rename?
}


ControlParameters::~ControlParameters() { }


void ControlParameters::set_mechanism(const std::string& mechanism_name)
{
    const Parameters* par = Parameters::get_parameters(mechanism_name);
    if (par == nullptr) return;
    this->par = par;
}


void ControlParameters::set_species(const std::vector<std::string> species_list, const std::vector<double> fractions_list)
{
    std::vector<index_t> species; species.reserve(species_list.size());
    for (const auto& species_name: species_list)
    {
        index_t index = par->get_species(species_name);
        if (index == par->invalid_index)
        {
            std::string message = "Invalid species (" + species_name + ") for mechanism " + par->mechanism_name;
            this->error_ID = LOG_ERROR(Error::severity::error, Error::type::preprocess, message, this->ID);
            return;
        }
        else if (index >= par->num_species)
        {
            std::string message = "Species index " + std::to_string(index) + " out of bound for mechanism " + par->mechanism_name;
            this->error_ID = LOG_ERROR(Error::severity::error, Error::type::preprocess, message, this->ID);
            return;
        }
        species.push_back(index);
    }

    this->set_species(species, fractions_list);
}


void ControlParameters::set_species(const std::vector<index_t>& species_list, const std::vector<double>& fractions_list)
{
    for (const auto& species: species_list)
        if (species == par->invalid_index)
        {
            std::string message = "Invalid species (" + std::to_string(species) + ") for mechanism " + par->mechanism_name;
            this->error_ID = LOG_ERROR(Error::severity::error, Error::type::preprocess, message, this->ID);
            return;
        }
        else if (species >= par->num_species)
        {
            std::string message = "Species index " + std::to_string(species) + " out of bound for mechanism " + par->mechanism_name;
            this->error_ID = LOG_ERROR(Error::severity::error, Error::type::preprocess, message, this->ID);
            return;
        }
    if (fractions_list.size() != species_list.size())
    {
        std::stringstream ss;
        ss << "The number of species and fractions must be equal: species_list=";
        ss << ::to_string((index_t*)species_list.data(), species_list.size());
        ss << ", fractions_list=";
        ss << ::to_string((double*)fractions_list.data(), fractions_list.size());
        this->error_ID = LOG_ERROR(Error::severity::error, Error::type::preprocess, ss.str(), this->ID);
        return;
    }
    if (fractions_list.size() > ControlParameters::max_species)
    {
        std::stringstream ss;
        ss << "Too many species: " << fractions_list.size() << " (" << ::to_string((index_t*)species_list.data(), species_list.size());
        ss << ", " << ::to_string((double*)fractions_list.data(), fractions_list.size()) << "). Perhaps try to change ControlParameters::max_species.";
        this->error_ID = LOG_ERROR(Error::severity::error, Error::type::preprocess, ss.str(), this->ID);
        return;
    }
    double sum_fraction = std::accumulate(fractions_list.begin(), fractions_list.end(), 0.0);
    if (std::abs(sum_fraction - 1.0) > 1.0e-10)
    {
        std::stringstream ss;
        ss << "The sum of fractions must be equal to 1.0, instead it is " << sum_fraction << ": fractions_list=";
        ss << ::to_string((double*)fractions_list.data(), fractions_list.size());
        ss << "]";
        this->error_ID = LOG_ERROR(Error::severity::error, Error::type::preprocess, ss.str(), this->ID);
        return;
    }

    this->num_initial_species = species_list.size();
    std::fill(this->species, this->species + ControlParameters::max_species, par->invalid_index);
    std::fill(this->fractions, this->fractions + ControlParameters::max_species, 0.0);
    std::copy(species_list.begin(), species_list.end(), species);
    std::copy(fractions_list.begin(), fractions_list.end(), fractions);
}


void ControlParameters::set_species(const std::initializer_list<std::string>& species_list, const std::initializer_list<double>& fractions_list)
{
    this->set_species(std::vector<std::string>(species_list), std::vector<double>(fractions_list));
}

void ControlParameters::set_species(const std::initializer_list<index_t>& species_list, const std::initializer_list<double>& fractions_list)
{
    this->set_species(std::vector<index_t>(species_list), std::vector<double>(fractions_list));
}


void ControlParameters::set_excitation_params(const std::initializer_list<double>& params_list)
{
    this->set_excitation_params(std::vector<double>(params_list));
}

void ControlParameters::set_excitation_params(const std::vector<double>& params_list)
{
    if (params_list.size() != Parameters::excitation_arg_nums[this->excitation_type])
    {
        std::stringstream ss;
        ss << "The number of excitation parameters must be equal to " << Parameters::excitation_arg_nums.at(this->excitation_type) << ": params_list=";
        ss << ::to_string((double*)params_list.data(), params_list.size()) << "; param_names={";
        ss << Parameters::excitation_arg_names.at(this->excitation_type) << "}; param_units={";
        ss << Parameters::excitation_arg_units.at(this->excitation_type) << "}";
        this->error_ID = LOG_ERROR(Error::severity::error, Error::type::preprocess, ss.str(), this->ID);
        return;
    }
    std::fill(this->excitation_params, this->excitation_params + ControlParameters::max_excitation_params, 0.0);
    std::copy(params_list.begin(), params_list.end(), excitation_params);
}


std::string ControlParameters::to_csv() const
{
    std::stringstream ss;
    auto format_double = [](std::ostream& os) -> std::ostream& {
        return os << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
    };
    if (par == nullptr)
    {
        return ",,,,,,,,,,,,,,,,,,,,,,,,";
    }

    ss << this->ID << "," << par->mechanism_name << "," << format_double << this->R_E << "," << format_double << this->ratio << ",";
    for (size_t index = 0; index < this->num_initial_species; ++index)
        ss << par->species_names[this->species[index]] << ";";
    ss << ",";
    for (size_t index = 0; index < this->num_initial_species; ++index)
        ss << format_double << this->fractions[index] << ";";
    ss << "," << format_double << this->P_amb << "," << format_double << this->T_inf << ",";
    ss << format_double << this->alpha_M << "," << format_double << this->P_v << ",";
    ss << format_double << this->mu_L << "," << format_double << this->rho_L << ",";
    ss << format_double << this->c_L << "," << format_double << this->surfactant << ",";
    ss << std::boolalpha << this->enable_heat_transfer << "," << std::boolalpha << this->enable_evaporation << ",";
    ss << std::boolalpha << this->enable_reactions << "," << std::boolalpha << this->enable_dissipated_energy << ",";
    ss << std::boolalpha << this->enable_van_der_waals << "," << std::boolalpha << this->enable_rate_thresholding << ",";
    ss << par->species_names[this->target_specie] << ",";
    ss << Parameters::excitation_names[this->excitation_type] << ",";
    for (size_t index = 0; index < Parameters::excitation_arg_nums[this->excitation_type]; ++index)
        ss << format_double << this->excitation_params[index] << ";";
    ss << "," << format_double << this->excitation_cycles << "," << format_double << this->ramp_up_cycles;

    return ss.str();
}


std::string ControlParameters::to_string(const bool with_code) const
{
    if (par == nullptr) return "Invalid mechanism name. ";

    std::stringstream ss;
    const size_t strw = 28;
    auto format_string = [](std::ostream& os) -> std::ostream& { return os << "    " << std::setw(strw); };
    auto format_bool = [](std::ostream& os) -> std::ostream& { return os << std::boolalpha; };
    auto format_double = [](std::ostream& os) -> std::ostream& {
        return os << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
    };
    auto species_to_string = [this, &with_code](const index_t species) -> std::string {
        if (species == this->par->invalid_index || species >= this->par->num_species)
        {
            if (with_code) return "\"X\"";
            else return "X";
        }

        if (with_code) return "\"" + this->par->species_names[species] + "\"";
        else return this->par->species_names[species];
    };
    std::vector<std::string> species_names;
    species_names.reserve(this->num_initial_species);
    for (size_t index = 0; index < this->num_initial_species; ++index)
        species_names.push_back(species_to_string(this->species[index]));

    if (with_code) ss << "ControlParameters::Builder{\n";
    ss << std::left;
    
    ss << format_string << ".ID"                         << " = " << this->ID << ",\n";
    ss << format_string << ".mechanism"                  << " = " << "\"" << par->mechanism_name << "\",\n";
    ss << format_string << ".R_E"                        << " = " << format_double << this->R_E                       << ",    // bubble equilibrium radius [m]\n";
    ss << format_string << ".ratio"                      << " = " << format_double << this->ratio                     << ",    // R_0/R_E for unforced oscillations [-]\n";
    ss << format_string << ".species"                    << " = " << ::to_string((std::string*)species_names.data(), this->num_initial_species)   << ",\n";
    ss << format_string << ".fractions"                  << " = " << ::to_string((double*)this->fractions, this->num_initial_species) << ",\n";
    ss << format_string << ".P_amb"                      << " = " << format_double << this->P_amb                     << ",    // ambient pressure [Pa]\n";
    ss << format_string << ".T_inf"                      << " = " << format_double << this->T_inf                     << ",    // ambient temperature [K]\n";
    ss << format_string << ".alpha_M"                    << " = " << format_double << this->alpha_M                   << ",    // water accommodation coefficient [-]\n";
    ss << format_string << ".P_v"                        << " = " << format_double << this->P_v                       << ",    // vapour pressure [Pa]\n";
    ss << format_string << ".mu_L"                       << " = " << format_double << this->mu_L                      << ",    // dynamic viscosity [Pa*s]\n";
    ss << format_string << ".rho_L"                      << " = " << format_double << this->rho_L                     << ",    // liquid density [kg/m^3]\n";
    ss << format_string << ".c_L"                        << " = " << format_double << this->c_L                       << ",    // sound speed [m/s]\n";
    ss << format_string << ".surfactant"                 << " = " << format_double << this->surfactant                << ",    // surface tension modifier [-]\n";
    ss << format_string << ".enable_heat_transfer"       << " = " << format_bool   << this->enable_heat_transfer      << ",\n";
    ss << format_string << ".enable_evaporation"         << " = " << format_bool   << this->enable_evaporation        << ",\n";
    ss << format_string << ".enable_reactions"           << " = " << format_bool   << this->enable_reactions          << ",\n";
    ss << format_string << ".enable_dissipated_energy"   << " = " << format_bool   << this->enable_dissipated_energy  << ",\n";
    ss << format_string << ".enable_van_der_waals"       << " = " << format_bool   << this->enable_van_der_waals      << ",\n";
    ss << format_string << ".enable_rate_thresholding"   << " = " << format_bool   << this->enable_rate_thresholding  << ",\n";
    ss << format_string << ".target_specie"              << " = " << species_to_string(this->target_specie)           << ",\n";
    ss << format_string << ".excitation_type"            << " = " << (with_code ? "Parameters::excitation::" : "") << Parameters::excitation_names[this->excitation_type] << "\n";
    ss << format_string << ".excitation_params"          << " = " << ::to_string((double*)this->excitation_params, Parameters::excitation_arg_nums[this->excitation_type]) << ",\n";
    ss << format_string << ".excitation_cycles"          << " = " << format_double << this->excitation_cycles         << ",    // number of excitation cycles [-]\n";
    ss << format_string << ".ramp_up_cycles"             << " = " << format_double << this->ramp_up_cycles            << "     // number of ramp-up cycles until full amplitude is reached [-]\n";

    if (with_code) ss << "}";
    ss << std::right;
    return ss.str();
}


ordered_json ControlParameters::to_json() const
{
    ordered_json j;
    if (par == nullptr) return j;
    std::vector<std::string> species_names;
    for (size_t index = 0; index < this->num_initial_species; ++index)
    {
        if (this->species[index] == par->invalid_index || this->species[index] >= par->num_species)
        {
            LOG_ERROR("Invalid species index: " + std::to_string(this->species[index]), this->ID);
            return j;
        }
        species_names.push_back(par->species_names.at(this->species[index]));
    }

    j["ID"] = this->ID;
    j["mechanism"] = this->par->mechanism_name;
    j["R_E"] = this->R_E;
    j["ratio"] = this->ratio;
    j["species"] = species_names;
    j["fractions"] = std::vector<double>(this->fractions, this->fractions + this->num_initial_species);
    j["P_amb"] = this->P_amb;
    j["T_inf"] = this->T_inf;
    j["alpha_M"] = this->alpha_M;
    j["P_v"] = this->P_v;
    j["mu_L"] = this->mu_L;
    j["rho_L"] = this->rho_L;
    j["c_L"] = this->c_L;
    j["surfactant"] = this->surfactant;
    j["enable_heat_transfer"] = this->enable_heat_transfer;
    j["enable_evaporation"] = this->enable_evaporation;
    j["enable_reactions"] = this->enable_reactions;
    j["enable_dissipated_energy"] = this->enable_dissipated_energy;
    j["enable_van_der_waals"] = this->enable_van_der_waals;
    j["enable_rate_thresholding"] = this->enable_rate_thresholding;
    j["target_specie"] = par->species_names.at(this->target_specie);
    j["excitation_type"] = Parameters::excitation_names.at(this->excitation_type);
    j["excitation_params"] = std::vector<double>(this->excitation_params, this->excitation_params + Parameters::excitation_arg_nums.at(this->excitation_type));
    j["excitation_cycles"] = this->excitation_cycles;
    j["ramp_up_cycles"] = this->ramp_up_cycles;

    return j;
}


std::ostream& operator<<(std::ostream& os, const ControlParameters& cpar)
{
    os << cpar.to_string(true);
    return os;
}


void ControlParameters::nondimensionalize(double &t, double* x) const
{
    t *= this->t_ref_inv;
    if (par == nullptr) return;
    if (x == nullptr) 
    {
        LOG_ERROR("x is nullptr");
        return;
    }

    const double R = x[0];
    const double V = 4.0 / 3.0 * std::numbers::pi * R * R * R;
    x[0] /= this->R_ref;
    x[1] *= this->t_ref / this->R_ref;
    x[2] /= this->T_ref;
    for (size_t i = 0; i < par->num_species; ++i)
    {
        const double n_i = x[i + 3] * V;
        const double n_dimless_i = n_i / this->c_ref;
        //const double n_log_i = std::log(n_dimless_i + this->epsilon);
        x[i+3] = n_log_i;
    }
    x[par->num_species + 3] /= this->E_diss_ref;
}


void ControlParameters::dimensionalize(double &t, double* x) const
{
    t *= this->t_ref;
    if (par == nullptr) return;
    if (x == nullptr) 
    {
        LOG_ERROR("x is nullptr");
        return;
    }

    const double R = x[0] * R_ref;
    const double V = 4.0 / 3.0 * std::numbers::pi * R * R * R;
    x[0] *= this->R_ref;
    x[1] *= this->R_ref * this->t_ref_inv;
    x[2] *= this->T_ref;
    for (size_t i = 0; i < par->num_species; ++i)
    {
        const double n_log_i = x[i + 3];
        const double n_dimless_i = std::exp(n_log_i) - this->epsilon;
        const double n_i = n_dimless_i * this->c_ref;
        const double c_i = n_i / V;
        x[i + 3] = c_i;
    }
    x[par->num_species + 3] *= this->E_diss_ref;
}


void ControlParameters::nondimensionalize_dot(double* x_dot, const double* x) const
{
    if (par == nullptr) return;
    if (x_dot == nullptr || x == nullptr) 
    {
        LOG_ERROR("x_dot or x is nullptr");
        return;
    }

    const double R = x[0];          // [m]
    const double R_dot = x[1];      // [m/s]
    const double V = 4.0 / 3.0 * std::numbers::pi * R * R * R;      // [m^3]
    const double V_dot = 4.0 * std::numbers::pi * R * R * R_dot;    // [m^3/s]
    x_dot[0] *= this->t_ref / this->R_ref;
    x_dot[1] *= this->t_ref * this->t_ref / this->R_ref;
    x_dot[2] *= this->t_ref / this->T_ref;
    for (index_t k = 0; k < par->num_species; ++k)
    {
        // x_dimless = ln(x * V / c_ref + epsilon) + log_offset = ln(f(x)) + log_offset
        // dx_dimless/dt_dimless = 1/f(x) * df/dt * dt/dt_dimless
            // f(x) = x * V / c_ref + epsilon
            // df/dt = (dx/dt * V + x * dV/dt) / c_ref
            // dt/dt_dimless = d/dt_dimless (t_dimless * t_ref) = t_ref
        // dx_dimless/dt_dimless = 1 / (x * V / c_ref + epsilon)  *  (dx/dt * V + x * dV/dt) / c_ref * t_ref
        const double fi = (x[k + 3] * V / this->c_ref + this->epsilon);
        const double dfdt = (x_dot[k + 3] * V + x[k + 3] * V_dot) / this->c_ref;
        x_dot[k + 3] = (1.0 / fi) * dfdt * this->t_ref;
    }
    x_dot[par->num_species + 3] *= this->t_ref / this->E_diss_ref;
}
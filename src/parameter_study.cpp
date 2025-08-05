#include <iomanip>
#include <sstream>
#include <fstream>
#include <numbers>
#include <cmath>
#include <regex>
#include <limits>
#include <filesystem>
#include <thread>

#include "nlohmann/json.hpp"
#include "parameter_study.h"
#include "ode_fun.h"
#include "ode_solver_sundials.h"

using ordered_json = nlohmann::ordered_json;


Range::Range():
    start(0.0),
    end(0.0),
    num_steps(0)
{ 
    LOG_ERROR("Range is not initialized witha any meaningful value.");
}


Range::Range(double start, double end, size_t num_steps):
    start(start),
    end(end),
    num_steps(num_steps)
{ }


std::string Range::to_array() const
{
    std::stringstream ss;
    ss << "{";

    size_t i_max;
    if (this->num_steps < 10)
        i_max = this->num_steps;
    else
        i_max = 7;
    for (size_t i = 0; i < i_max; ++i)
    {
        ss << this->operator[](i);
        if (i < i_max - 1)
        {
            ss << ", ";
        }
    }
    if (this->num_steps >= 10)
    {
        ss << ", ..., ";
        ss << this->operator[](this->num_steps - 2);
        ss << ", ";
        ss << this->operator[](this->num_steps - 1);
    }
    ss << "}";
    return ss.str();
}


std::ostream &operator<<(std::ostream &os, const Range &range)
{
    os << range.to_string();
    return os;
}


size_t Range::get_num_steps() const
{
    return this->num_steps;
}


Const::Const(double value):
    Range(value, value, 1)
{}


double Const::operator[](size_t i) const
{
    (void)i;
    return this->start;
}


std::string Const::to_string() const
{
    return "Const(" + std::to_string(this->start) + ")";
}


ordered_json Const::to_json() const
{
    return {
        {"type", "Const"},
        {"value", this->start}
    };
}


LinearRange::LinearRange(double start, double end, size_t num_steps):
    Range(start, end, std::max(num_steps, (size_t)1))
{}


double LinearRange::operator[](size_t i) const
{
    if (this->num_steps <= 1)
    {
        return this->start;
    }
    if (i >= this->num_steps)
    {
        return this->end;
    }

    return this->start + i * (this->end - this->start) / (this->num_steps - 1);
}


std::string LinearRange::to_string() const
{
    return "LinearRange(" + std::to_string(this->start) + ", " + std::to_string(this->end) + ", " + std::to_string(this->num_steps) + ")";
}


ordered_json LinearRange::to_json() const
{
    return {
        {"type", "LinearRange"},
        {"start", this->start},
        {"end", this->end},
        {"num_steps", this->num_steps}
    };
}


PowRange::PowRange(double start, double end, size_t num_steps, double base):
    Range(start, end, std::max(num_steps, (size_t)1))
{
    if (base <= 0.0)
    {
        LOG_ERROR("The value of base (" + std::to_string(base) + ") must be positive. Otherwise the formula f(x) = a * x^base + c breaks.");
        this->start = 1.0;
    }
    this->base = base;
    this->a = (end - start) / std::pow(num_steps - 1, base);
    this->c = start;
}


double PowRange::operator[](size_t i) const
{
    if (this->num_steps <= 1)
    {
        return this->start;
    }
    if (i >= this->num_steps)
    {
        return this->end;
    }

    return this->a * std::pow(i, this->base) + this->c;
}


std::string PowRange::to_string() const
{
    return "PowRange(" + std::to_string(this->start) + ", " + std::to_string(this->end) + ", " + std::to_string(this->num_steps) + ", " + std::to_string(this->base) + ")";
}


ordered_json PowRange::to_json() const
{
    return {
        {"type", "PowRange"},
        {"start", this->start},
        {"end", this->end},
        {"num_steps", this->num_steps},
        {"base", this->base}
    };
}


std::unique_ptr<Range> get_unique_ptr(const ParameterCombinator::AnyRange range)
{
    if (std::holds_alternative<Const>(range))
    {
        return std::make_unique<Const>(std::get<Const>(range));
    }
    else if (std::holds_alternative<LinearRange>(range))
    {
        return std::make_unique<LinearRange>(std::get<LinearRange>(range));
    }
    else if (std::holds_alternative<PowRange>(range))
    {
        return std::make_unique<PowRange>(std::get<PowRange>(range));
    }
    return nullptr;
}


void ParameterCombinator::init(const ParameterCombinator::Builder &builder)
{
    this->combination_ID =              0;
    this->total_combination_count =     1;
    this->mechanism =                   builder.mechanism;
    this->R_E =                         get_unique_ptr(builder.R_E);
    this->species =                     builder.species;
    this->fractions =                   builder.fractions;
    this->P_amb =                       get_unique_ptr(builder.P_amb);
    this->T_inf =                       get_unique_ptr(builder.T_inf);
    this->alfa_M =                      get_unique_ptr(builder.alfa_M);
    this->P_v =                         get_unique_ptr(builder.P_v);
    this->mu_L =                        get_unique_ptr(builder.mu_L);
    this->rho_L =                       get_unique_ptr(builder.rho_L);
    this->c_L =                         get_unique_ptr(builder.c_L);
    this->surfactant =                  get_unique_ptr(builder.surfactant);
    this->enable_heat_transfer =        builder.enable_heat_transfer;
    this->enable_evaporation =          builder.enable_evaporation;
    this->enable_reactions =            builder.enable_reactions;
    this->enable_dissipated_energy =    builder.enable_dissipated_energy;
    this->target_specie =               builder.target_specie;
    this->excitation_type =             builder.excitation_type;

    this->excitation_params.reserve(builder.excitation_params.size());
    for (const auto &excitation_param : builder.excitation_params)
    {
        this->excitation_params.push_back(get_unique_ptr(excitation_param));
    }

    this->total_combination_count *= this->R_E->get_num_steps();
    this->total_combination_count *= this->P_amb->get_num_steps();
    this->total_combination_count *= this->T_inf->get_num_steps();
    this->total_combination_count *= this->alfa_M->get_num_steps();
    this->total_combination_count *= this->P_v->get_num_steps();
    this->total_combination_count *= this->mu_L->get_num_steps();
    this->total_combination_count *= this->rho_L->get_num_steps();
    this->total_combination_count *= this->c_L->get_num_steps();
    this->total_combination_count *= this->surfactant->get_num_steps();
    for (const auto &excitation_param : this->excitation_params)
    {
        this->total_combination_count *= excitation_param->get_num_steps();
    }
}


ParameterCombinator::ParameterCombinator(const ParameterCombinator::Builder &builder)
{
    this->init(builder);
}


ParameterCombinator::AnyRange get_range(const ordered_json& j, ParameterCombinator::AnyRange default_range)
{
    if (!j.is_object())
    {
        LOG_ERROR(
            "Expected JSON object, instead found " + j.dump() + \
            ". Available options: Const(value), LinearRange(start, end, num_steps), PowRange(start, end, num_steps, base). " + \
            "E.g.: {\"type\": \"LinearRange\", \"start\": 0.0, \"end\": 1.0, \"num_steps\": 10}"
        );
        return default_range;
    }
    if (!j.contains("type") || !j.at("type").is_string())
    {
        LOG_ERROR(
            "Invalid \"type\" in JSON object " + j.dump() + \
            ". Available options: Const(value), LinearRange(start, end, num_steps), PowRange(start, end, num_steps, base). " + \
            "e.g.: {\"type\": \"LinearRange\", \"start\": 0.0, \"end\": 1.0, \"num_steps\": 10}"
        );
        return default_range;
    }
    std::string type = j.at("type").get<std::string>();
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    std::replace(type.begin(), type.end(), '_', ' ');
    type.erase(std::remove_if(type.begin(), type.end(), [](unsigned char c) { return std::isspace(c); }), type.end());
    
    if (type == "const")
    {
        if (!j.contains("value") || !j.at("value").is_number())
        {
            LOG_ERROR(
                "Invalid arguments for JSON object " + j.dump() + \
                ". For Const, you must provide a value of type double. " + \
                "e.g.: {\"type\": \"Const\", \"value\": 0.0}"
            );
            return default_range;
        }
        return Const(j.at("value").get<double>());
    }
    else if (type == "linearrange")
    {
        if (
            !j.contains("start") || !j.contains("end") || !j.contains("num_steps")
            || !j.at("start").is_number() || !j.at("end").is_number() || !j.at("num_steps").is_number()
        )
        {
            LOG_ERROR(
                "Invalid arguments for JSON object " + j.dump() + \
                ". For LinearRange, you must provide start, end of type double and and num_steps of type integer. " + \
                "e.g.: {\"type\": \"LinearRange\", \"start\": 0.0, \"end\": 1.0, \"num_steps\": 10}"
            );
            return default_range;
        }
        return LinearRange(
            j.at("start").get<double>(),
            j.at("end").get<double>(),
            j.at("num_steps").get<size_t>()
        );
    }
    else if (type == "powrange")
    {
        if (
            !j.contains("start") || !j.contains("end") || !j.contains("num_steps")
            || !j.at("start").is_number() || !j.at("end").is_number() || !j.at("num_steps").is_number()
            || (j.contains("base") && !j.at("base").is_number())
        )
        {
            LOG_ERROR(
                "Invalid arguments for JSON onject " + j.dump() + \
                ". For PowRange, you must provide start, end, and optionally base of type double. As well as num_steps of type integer. " + \
                "e.g.: {\"type\": \"PowRange\", \"start\": 0.0, \"end\": 1.0, \"num_steps\": 10, \"base\": 2.0}"
            );
            return default_range;
        }
        return PowRange(
            j.at("start").get<double>(),
            j.at("end").get<double>(),
            j.at("num_steps").get<size_t>(),
            j.value("base", 2.0)
        );
    }
    else
    {
        LOG_ERROR("Invalid Range type: " + type + \
            ". Available options: Const(value), LinearRange(start, end, num_steps), PowRange(start, end, num_steps, base). " + \
            "e.g.: {\"type\": \"LinearRange\", \"start\": 0.0, \"end\": 1.0, \"num_steps\": 10}"
        );
        return default_range;
    }
}


ParameterCombinator::AnyRange get_range(const ordered_json& j, const std::string &key, ParameterCombinator::AnyRange default_range)
{
    if (!j.contains(key))
    {
        LOG_ERROR(
            Error::severity::warning,
            Error::type::preprocess,
            "Key \"" + key + "\" not found in JSON object. " + \
            "Using default value. "
        );
        return default_range;
    }
    else
    {
        return get_range(j.at(key), default_range);
    }
}


template<typename T>
T get_value(const ordered_json& j, const std::string& key, const T& default_value)
{
    if constexpr (std::is_same_v<T, std::vector<double>> || std::is_same_v<T, std::vector<std::string>>)
    {
        if (j.contains(key) && !j.at(key).is_array())
        {
            LOG_ERROR(
                Error::severity::warning,
                Error::type::preprocess,
                "Expected JSON array for key \"" + key + "\", instead found " + j.at(key).dump() + ". Using default value."
            );
            return default_value;
        }
    }
    
    if (j.contains(key))
    {
        return j.at(key).get<T>();
    }

    LOG_ERROR(
        Error::severity::warning,
        Error::type::preprocess,
        "Key \"" + key + "\" not found in JSON object. Using default value. "
    );
    return default_value;
}


void ParameterCombinator::init(const ordered_json& j)
{
    auto builder = ParameterCombinator::Builder{};
    try
    {
        builder.mechanism = Parameters::string_to_mechanism(
                                            get_value<std::string>              (j, "mechanism",                Parameters::mechanism_names.at(builder.mechanism))
        );
        builder.R_E =                       get_range(j, "R_E",                      builder.R_E);
        builder.species =                   get_value<std::vector<std::string>> (j, "species",                  builder.species);
        builder.fractions =                 get_value<std::vector<double>>      (j, "fractions",                builder.fractions);
        builder.P_amb =                     get_range                           (j, "P_amb",                    builder.P_amb);
        builder.T_inf =                     get_range                           (j, "T_inf",                    builder.T_inf);
        builder.alfa_M =                    get_range                           (j, "alfa_M",                   builder.alfa_M);
        builder.P_v =                       get_range                           (j, "P_v",                      builder.P_v);
        builder.mu_L =                      get_range                           (j, "mu_L",                     builder.mu_L);
        builder.rho_L =                     get_range                           (j, "rho_L",                    builder.rho_L);
        builder.c_L =                       get_range                           (j, "c_L",                      builder.c_L);
        builder.surfactant =                get_range                           (j, "surfactant",               builder.surfactant);
        builder.enable_heat_transfer =      get_value<bool>                     (j, "enable_heat_transfer",     builder.enable_heat_transfer);
        builder.enable_evaporation =        get_value<bool>                     (j, "enable_evaporation",       builder.enable_evaporation);
        builder.enable_reactions =          get_value<bool>                     (j, "enable_reactions",         builder.enable_reactions);
        builder.enable_dissipated_energy =  get_value<bool>                     (j, "enable_dissipated_energy", builder.enable_dissipated_energy);
        builder.target_specie =             get_value<std::string>              (j, "target_specie",            builder.target_specie);
        builder.excitation_type = Parameters::string_to_excitation(
                                            get_value<std::string>              (j, "excitation_type",          Parameters::excitation_names.at(builder.excitation_type))
        );

        if (j.contains("excitation_params"))
        {
            if (!j.at("excitation_params").is_array())
            {
                LOG_ERROR(
                    "\"excitation_params\" should be a list, instead it is " + j.at("excitation_params").dump() + \
                    ". E.g.: {\"excitation_params\": [{\"type\": \"Const\", \"value\": 0.0}, ...]}"
                );
            }
            else if (j.at("excitation_params").size() != Parameters::excitation_arg_nums.at(builder.excitation_type))
            {
                LOG_ERROR(
                    "Invalid number of excitation parameters for excitation type " + std::string(Parameters::excitation_names.at(builder.excitation_type)) + \
                    ". Expected " + std::to_string(Parameters::excitation_arg_nums.at(builder.excitation_type)) + \
                    ", but got " + std::to_string(j.at("excitation_params").size()) + \
                    ". Excitation params are: " + Parameters::excitation_arg_names.at(builder.excitation_type)
                );
            }
            else
            {
                builder.excitation_params.resize(j.at("excitation_params").size(), Const(0.0));
                for (size_t i=0; i < j.at("excitation_params").size(); ++i)
                {
                    ordered_json param = j.at("excitation_params").at(i);
                    if (param != nullptr)
                        builder.excitation_params.at(i) = get_range(param, builder.excitation_params.at(i));
                    
                }
            }
        }
    }
    catch(const std::exception& e)
    {
        LOG_ERROR(
            "Error parsing JSON file: " + std::string(e.what())
        );
    }

    this->init(builder);
}


ParameterCombinator::ParameterCombinator(const ordered_json& j)
{
    this->init(j);
}



ParameterCombinator::ParameterCombinator(const std::string& json_path)
{
    // Open JSON file
    std::ifstream input_file(json_path);
    if (!input_file.is_open())
    {
        this->init(ParameterCombinator::Builder());
        LOG_ERROR(
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
        this->init(ParameterCombinator::Builder());
        LOG_ERROR(
            "Error parsing JSON file: " + std::string(e.what())
        );
        return;
    }
    input_file.close();

    // Parse JSON data
    if (!j.contains("parameter_study"))
    {
        this->init(ParameterCombinator::Builder());
        LOG_ERROR(
            "JSON file does not contain 'parameter_study' key: " + j.dump()
        );
        return;
    }

    this->init(j.at("parameter_study"));
}


std::string ParameterCombinator::to_string(const bool with_code) const
{
    const Parameters* par = Parameters::get_parameters(this->mechanism);
    if (par == nullptr)
    {
        return "";
    }
    std::stringstream ss;
    const size_t strw = 28;
    const size_t rangew = 48;
    auto format_field_name = [](std::ostream& os) -> std::ostream& {
        return os << "    " << std::setw(strw);
    };
    auto format_range = [](Range *range, const bool last = false) -> std::string {
        std::stringstream ss;
        ss << std::left << std::setw(rangew) << (range->to_string() + (last ? " " : ",")) << "    // " << range->to_array();
        return ss.str();
    };
    auto species_to_string = [&with_code](const std::string name) -> std::string {
        if (with_code) return "\"" + name + "\"";
        else return name;
    };
    std::vector<std::string> species_names;
    species_names.reserve(this->species.size());
    for (auto &species_name : this->species)
    {
        species_names.push_back(species_to_string(species_name));
    }

    if (with_code) ss << "ParameterCombinator::Builder{\n";
    ss << std::left;

    ss << format_field_name << ".mechanism" << " = " << ((with_code ? "Parameters::mechanism::" : "") + par->model) << ",\n";
    ss << format_field_name << ".R_E" << " = " << format_range(this->R_E.get()) << "\n";
    ss << format_field_name << ".species" << " = " << ::to_string((std::string*)species_names.data(), species_names.size()) << ",\n";
    ss << format_field_name << ".fractions" << " = " << ::to_string((double*)this->fractions.data(), this->fractions.size()) << ",\n";
    ss << format_field_name << ".P_amb" << " = " << format_range(this->P_amb.get()) << "\n";
    ss << format_field_name << ".T_inf" << " = " << format_range(this->T_inf.get()) << "\n";
    ss << format_field_name << ".alfa_M" << " = " << format_range(this->alfa_M.get()) << "\n";
    ss << format_field_name << ".P_v" << " = " << format_range(this->P_v.get()) << "\n";
    ss << format_field_name << ".mu_L" << " = " << format_range(this->mu_L.get()) << "\n";
    ss << format_field_name << ".rho_L" << " = " << format_range(this->rho_L.get()) << "\n";
    ss << format_field_name << ".c_L" << " = " << format_range(this->c_L.get()) << "\n";
    ss << format_field_name << ".surfactant" << " = " << format_range(this->surfactant.get()) << "\n";
    ss << format_field_name << ".enable_heat_transfer" << " = " << std::boolalpha << this->enable_heat_transfer << ",\n";
    ss << format_field_name << ".enable_evaporation" << " = " << std::boolalpha << this->enable_evaporation << ",\n";
    ss << format_field_name << ".enable_reactions" << " = " << std::boolalpha << this->enable_reactions << ",\n";
    ss << format_field_name << ".enable_dissipated_energy" << " = " << std::boolalpha << this->enable_dissipated_energy << ",\n";
    ss << format_field_name << ".target_specie" << " = " << species_to_string(this->target_specie) << ",\n";
    ss << format_field_name << ".excitation_params" << " = {";
    for (size_t i = 0; i < this->excitation_params.size(); ++i)
    {
        if (i == 0)
            ss << "\n";
        ss << "        " << format_range(this->excitation_params.at(i).get(), i == this->excitation_params.size() - 1);
        ss << "\n";
    }

    ss << "    },\n";
    ss << format_field_name << ".excitation_type"            << " = " << (with_code ? "Parameters::excitation::" : "") << Parameters::excitation_names[this->excitation_type] << "\n";

    if (with_code) ss << "}";
    ss << std::right;
    return ss.str();
}


ordered_json ParameterCombinator::to_json() const
{
    ordered_json j;

    j["mechanism"] = Parameters::mechanism_names.at(this->mechanism);
    j["R_E"] = this->R_E->to_json();
    j["species"] = this->species;
    j["fractions"] = this->fractions;
    j["P_amb"] = this->P_amb->to_json();
    j["T_inf"] = this->T_inf->to_json();
    j["alfa_M"] = this->alfa_M->to_json();
    j["P_v"] = this->P_v->to_json();
    j["mu_L"] = this->mu_L->to_json();
    j["rho_L"] = this->rho_L->to_json();
    j["c_L"] = this->c_L->to_json();
    j["surfactant"] = this->surfactant->to_json();
    j["enable_heat_transfer"] = this->enable_heat_transfer;
    j["enable_evaporation"] = this->enable_evaporation;
    j["enable_reactions"] = this->enable_reactions;
    j["enable_dissipated_energy"] = this->enable_dissipated_energy;
    j["target_specie"] = this->target_specie;
    j["excitation_params"] = ordered_json::array();
    for (const auto& excitation_param : this->excitation_params)
    {
        j["excitation_params"].push_back(excitation_param->to_json());
    }
    j["excitation_type"] = Parameters::excitation_names.at(this->excitation_type);

    return {{"parameter_study", j}};
}


std::ostream &operator<<(std::ostream &os, const ParameterCombinator &pc)
{
    os << pc.to_string(true);
    return os;
}


size_t ParameterCombinator::get_total_combination_count() const
{
    return this->total_combination_count;
}


size_t ParameterCombinator::get_next_combination_ID() const
{
    return this->combination_ID.load();
}


std::pair<is_success, ControlParameters> ParameterCombinator::get_next_combination()
{
    const size_t combination_ID = this->combination_ID.fetch_add(1, std::memory_order_relaxed);    
    size_t prod = 1;

    auto get_value = [&combination_ID, &prod](const std::unique_ptr<Range> &range) -> double {
        double value = range->operator[](combination_ID / prod % range->get_num_steps());
        prod *= range->get_num_steps();
        return value;
    };

    std::vector<double> excitation_values(this->excitation_params.size());
    for (size_t i = 0; i < this->excitation_params.size(); ++i)
    {
        excitation_values.at(i) = get_value(this->excitation_params.at(i));
    }

    ControlParameters cpar{ControlParameters::Builder{
        .ID = combination_ID,
        .mechanism = this->mechanism,
        .R_E = get_value(this->R_E),
        .species = this->species,
        .fractions = this->fractions,
        .P_amb = get_value(this->P_amb),
        .T_inf = get_value(this->T_inf),
        .alfa_M = get_value(this->alfa_M),
        .P_v = get_value(this->P_v),
        .mu_L = get_value(this->mu_L),
        .rho_L = get_value(this->rho_L),
        .c_L = get_value(this->c_L),
        .surfactant = get_value(this->surfactant),
        .enable_heat_transfer = this->enable_heat_transfer,
        .enable_evaporation = this->enable_evaporation,
        .enable_reactions = this->enable_reactions,
        .enable_dissipated_energy = this->enable_dissipated_energy,
        .target_specie = this->target_specie,
        .excitation_params = excitation_values,
        .excitation_type = this->excitation_type
    }};

    is_success success = combination_ID < this->total_combination_count;

    return {success, cpar};
}


const Parameters* ParameterCombinator::get_mechanism_parameters() const
{
    const Parameters* par = Parameters::get_parameters(this->mechanism);
    if (par == nullptr)
    {
        LOG_ERROR("Mechanism " + std::to_string(this->mechanism) + " is not found.");
    }
    return par;
}


double get_dissipated_energy(const OdeSolution &sol)
{
    if (sol.x.empty()) return 0.0;
    if (sol.x.back().empty()) return 0.0;
    return sol.x.back().back();
}

double get_n_target(const OdeSolution &sol, const ControlParameters &cpar)
{
    if (sol.x.empty()) return 0.0;
    if (sol.x.back().empty()) return 0.0;
    const Parameters* par = Parameters::get_parameters(cpar.mechanism);
    if (par == nullptr) return 0.0;
    if (cpar.target_specie == par->invalid_index) return 0.0;
    if (sol.x.front().size() != size_t(4 + par->num_species)) return 0.0;
    if (sol.x.back().size() != size_t(4 + par->num_species)) return 0.0;

    const double R_last = 100.0 * sol.x.back()[0];  // [cm]
    const double V_last = 4.0 / 3.0 * std::numbers::pi * std::pow(R_last, 3); // [cm^3]
    const double c_target = sol.x.back()[3+cpar.target_specie];  // [mol/cm^3]

    return c_target * V_last;  // [mol]
}

double get_energy_demand(const OdeSolution &sol, const ControlParameters &cpar)
{
    const Parameters* par = Parameters::get_parameters(cpar.mechanism);
    if (par == nullptr) return SimulationData::infinite_energy_demand;
    if (cpar.target_specie == par->invalid_index) return SimulationData::infinite_energy_demand;

    const double dissipated_energy = get_dissipated_energy(sol);    // [J]
    const double n_target = get_n_target(sol, cpar);    // [mol]
    const double m_target = 1.0e-3 * n_target * par->W[cpar.target_specie];  // [kg]

    if (dissipated_energy < -1.0e-8)
        LOG_ERROR(Error::severity::warning, Error::type::postprocess, "Dissipated energy is negative: " + std::to_string(dissipated_energy), cpar.ID);
    if (n_target < -1.0e-8)
        LOG_ERROR(Error::severity::warning, Error::type::postprocess, "Target specie concentration is negative: " + std::to_string(n_target), cpar.ID);

    if (m_target < 10*std::numeric_limits<double>::min())
        return SimulationData::infinite_energy_demand;

    double energy_demand = 1.0e-6 * dissipated_energy / m_target;  // [MJ/kg]
    if (energy_demand < 0.0)
    {
        LOG_ERROR(Error::severity::warning, Error::type::postprocess, "Energy demand is negative: " + std::to_string(energy_demand), cpar.ID);
        return SimulationData::infinite_energy_demand;
    }

    return energy_demand;
}


const std::string SimulationData::csv_header = std::string("dissipated_energy,n_target_specie,energy_demand,")
                                             + std::string(ControlParameters::csv_header) + std::string(",")
                                             + std::string(OdeSolution::csv_header) + std::string(",") + std::string(Error::csv_header);

const Error SimulationData::no_error = Error(Error::severity::info, Error::type::general, "No error", "", __FILE__, __LINE__, 0);

const double SimulationData::infinite_energy_demand = std::numeric_limits<double>::infinity();

SimulationData::SimulationData(const ControlParameters &cpar, const OdeSolution &sol):
    cpar(cpar),    
    sol(sol),
    dissipated_energy(get_dissipated_energy(sol)),
    n_target_specie(get_n_target(sol, cpar)),
    energy_demand(get_energy_demand(sol, cpar))
{}


std::string SimulationData::to_csv() const
{
    std::stringstream ss;
    auto format_double = [](std::ostream& os) -> std::ostream& {
        return os << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
    };

    ss << format_double << this->dissipated_energy << ",";
    ss << format_double << this->n_target_specie << ",";
    ss << format_double << this->energy_demand << ",";
    ss << this->cpar.to_csv() << ",";
    ss << this->sol.to_csv() << ",";

    if (this->sol.error_ID == ErrorHandler::no_error)
    {
        ss << this->no_error.to_csv();
    } else
    {
        Error error = ErrorHandler::get_error(this->sol.error_ID);
        ss << error.to_csv();
    }


    return ss.str();
}


std::string SimulationData::to_string() const
{
    const size_t strw = 28;
    auto format_double = [](std::ostream& os) -> std::ostream& {
        return os << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
    };

    std::stringstream ss;
    ss << std::left;
    ss << "SimulationData{\n";
    ss << std::setw(strw) << "    .dissipated_energy"  << " = " << format_double << this->dissipated_energy    << ",    // [J]\n";
    ss << std::setw(strw) << "    .n_target_specie"    << " = " << format_double << this->n_target_specie      << ",    // [mol]\n";
    ss << std::setw(strw) << "    .energy_demand"      << " = " << format_double << this->energy_demand        << ",    // [MJ/kg]\n";
    ss << "    .cpar = ControlParameters{";
    ss << std::regex_replace(this->cpar.to_string(true), std::regex("\n"), "\n    ");
    ss << "},\n    .sol = ";
    ss << std::regex_replace(this->sol.to_string(false, true), std::regex("\n"), "\n    ");
    ss << ",\n}" << std::right;

    return ss.str();
}


std::string SimulationData::to_small_string(const ParameterCombinator &ps, const double best_energy_demand, const bool colored) const
{
    std::stringstream ss;
    std::string total_combinations = std::to_string(ps.get_total_combination_count());
    auto format_double = [](std::ostream& os) -> std::ostream& {
        return os << std::scientific << std::setprecision(4);
    };

    ss << std::setw(total_combinations.size()) << std::to_string(this->cpar.ID) << "/" << total_combinations << ": ";
    if (colored)
        ss << colors::bold << (this->sol.success() ? colors::green : colors::red) << (this->sol.success() ? "success" : " failed") << colors::reset << "; ";
    else
        ss << (this->sol.success() ? "success" : " failed") << "; ";
    ss << "runtime=" << std::setw(10) << Timer::format_time(this->sol.runtime) << "; ";
    ss << "num_steps=" << std::setw(8) << this->sol.num_steps << ";   |   ";

    if (ps.R_E->get_num_steps() > 1)
        ss << "R_E=" << format_double << this->cpar.R_E << " m; ";
    if (ps.P_amb->get_num_steps() > 1)
        ss << "P_amb=" << format_double << this->cpar.P_amb << " Pa; ";
    if (ps.T_inf->get_num_steps() > 1)
        ss << "T_inf=" << format_double << this->cpar.T_inf << " K; ";
    if (ps.alfa_M->get_num_steps() > 1)
        ss << "alfa_M=" << format_double << this->cpar.alfa_M << "; ";
    if (ps.P_v->get_num_steps() > 1)
        ss << "P_v=" << format_double << this->cpar.P_v << " Pa; ";
    if (ps.mu_L->get_num_steps() > 1)
        ss << "mu_L=" << format_double << this->cpar.mu_L << " Pa*s; ";
    if (ps.rho_L->get_num_steps() > 1)
        ss << "rho_L=" << format_double << this->cpar.rho_L << " kg/m^3; ";
    if (ps.c_L->get_num_steps() > 1)
        ss << "c_L=" << format_double << this->cpar.c_L << " m/s; ";
    if (ps.surfactant->get_num_steps() > 1)
        ss << "surfactant=" << format_double << this->cpar.surfactant << "; ";

    std::string excitation_arg_names = Parameters::excitation_arg_names.at(ps.excitation_type);
    std::string excitation_arg_units = Parameters::excitation_arg_units.at(ps.excitation_type);
    auto name_begin = excitation_arg_names.begin();
    auto unit_begin = excitation_arg_units.begin();
    auto name_end = excitation_arg_names.end();
    auto unit_end = excitation_arg_units.end();
    for (size_t i = 0; i < ps.excitation_params.size(); ++i)
    {
        name_end = std::find(name_begin, excitation_arg_names.end(), ' ');
        unit_end = std::find(unit_begin, excitation_arg_units.end(), ' ');
        std::string excitation_arg_name = std::string(name_begin, name_end);
        std::string excitation_arg_unit = std::string(unit_begin, unit_end);
        if (ps.excitation_params[i]->get_num_steps() > 1)
            ss << excitation_arg_name << "=" << format_double << this->cpar.excitation_params[i] << " " << excitation_arg_unit << "; ";
        name_begin = name_end + 1;
        unit_begin = unit_end + 1;
    }

    ss << "   |   ";
    ss << "energy_demand=" << format_double << this->energy_demand << " MJ/kg ";
    if (colored)
        ss << colors::bold << "(best=" << format_double << best_energy_demand << " MJ/kg)" << colors::reset;
    else
        ss << "(best=" << format_double << best_energy_demand << " MJ/kg)";

    return ss.str();
}


// Helper function to split a space-separated string into a vector of strings
std::vector<std::string> split_string(const std::string& str) {
    std::vector<std::string> result;
    std::stringstream sstream(str);
    std::string token;
    while (sstream >> token) {
        result.push_back(token);
    }
    return result;
}


ordered_json SimulationData::to_json() const
{
    ordered_json j;
    j["dissipated_energy"] = this->dissipated_energy;
    j["n_target_specie"] = this->n_target_specie;
    j["energy_demand"] = this->energy_demand;
    j["cpar"] = this->cpar.to_json();
    j["sol"] = this->sol.to_json();
    j["excitation"] = ordered_json::object({
        {"type", Parameters::excitation_names.at(this->cpar.excitation_type)},
        {
            "names",
            split_string(
                Parameters::excitation_arg_names.at(this->cpar.excitation_type)
            )
        },
        {
            "units",
            split_string(
                Parameters::excitation_arg_units.at(this->cpar.excitation_type)
            )
        }
    });

    const Parameters* par = Parameters::get_parameters(this->cpar.mechanism);
    if (par == nullptr)
    {
        LOG_ERROR("Mechanism " + std::to_string(this->cpar.mechanism) + " is not found.");
        return j;
    }
    j["mechanism"] = ordered_json::object({
        {"model", par->model},
        {"num_species", par->num_species},
        {"num_reactions", par->num_reactions},
        {"species_names", par->species_names}
    });
    j["version"] = VERSION;

    return j;
}


void SimulationData::save_json_with_binary(const std::string &json_path) const
{
    // Open file as text
    std::ofstream output_file(json_path, std::ios::out);
    if (!output_file.is_open())
    {
        LOG_ERROR("Could not open output JSON file: " + json_path);
        return;
    }

    // Save JSON data + <BINARY> marker
    ordered_json j = this->to_json();
    output_file << std::setw(4) << j << std::endl;
    output_file << std::endl << "<BINARY>";
    output_file.close();

    // Open file as binary
    std::ofstream binary_output_file(json_path, std::ios::app | std::ios::binary);
    if (!binary_output_file.is_open())
    {
        LOG_ERROR("Could not open output file as binary: " + json_path);
        return;
    }

    // Save sol.t (1D array)
    if (this->sol.t.size() != this->sol.x.size())
    {
        LOG_ERROR("Mismatch between sol.t.size() and sol.x.size(): " + std::to_string(this->sol.t.size()) + " != " + std::to_string(this->sol.x.size()));
        return;
    }
    binary_output_file.write(reinterpret_cast<const char*>(this->sol.t.data()), this->sol.t.size() * sizeof(double));

    // Save sol.x (2D array)
    for (const auto& row : this->sol.x)
    {
        if (row.size() != this->sol.num_dim)
        {
            LOG_ERROR("Mismatch between sol.x[].size() and sol.num_dim: " + std::to_string(row.size()) + " != " + std::to_string(this->sol.num_dim));
            return;
        }
        binary_output_file.write(reinterpret_cast<const char*>(row.data()), row.size() * sizeof(double));
    }
    binary_output_file.close();
}


std::ostream &operator<<(std::ostream &os, const SimulationData &data)
{
    os << data.to_string();
    return os;
}


// Verify, that the path provided as string is valid. Automatically count folders with the same base name. 
// Create a new folder, return the path to the new folder.
std::string verify_save_folder(const std::string save_folder)
{
    std::filesystem::path save_folder_path(save_folder);

    // count existing folders with the same base name
    size_t folder_count = 0;
    if (std::filesystem::exists(save_folder_path.parent_path()))
    {
        for (const auto& entry : std::filesystem::directory_iterator(save_folder_path.parent_path()))
        {
            if (entry.is_directory() && entry.path().filename().string().find(save_folder_path.filename().string()) == 0)
            {
                folder_count++;
            }
        }
    }
    
    // create folder name
    save_folder_path = save_folder_path.parent_path() / (save_folder_path.filename().string() + std::to_string(folder_count + 1));
    if (std::filesystem::exists(save_folder_path))
    {
        LOG_ERROR("Save folder already exists: " + save_folder_path.string());
        return "";
    }

    // create save folder
    if (!std::filesystem::create_directories(save_folder_path))
    {
        LOG_ERROR("Failed to create save folder: " + save_folder_path.string());
        return "";
    }

    return save_folder_path.string();;
}


ParameterStudy::ParameterStudy(
    ParameterCombinator &parameter_combinator,
    std::string save_folder,
    const double t_max,
    const double timeout
):
    parameter_combinator(parameter_combinator),
    save_folder(verify_save_folder(save_folder)),
    best_energy_demand(SimulationData::infinite_energy_demand),
    t_max(t_max),
    timeout(timeout),
    successful_simulations(0),
    total_simulations(0)
{
    // save general information
    if (this->save_folder.empty()) return;
    std::filesystem::path save_folder_path(this->save_folder);
    std::filesystem::path general_info_file_path = save_folder_path / "bruteforce_parameter_study_settings.txt";
    std::ofstream general_info_file(general_info_file_path);
    if (!general_info_file.is_open())
    {
        LOG_ERROR("Failed to open file: " + general_info_file_path.string());
    } else {
        const Parameters* par = Parameters::get_parameters(parameter_combinator.mechanism);
        general_info_file << "version: " << VERSION << "\n";
        general_info_file << "datetime: " << Timer::current_time() << "\n";
        general_info_file << "total_combination_count: " << parameter_combinator.get_total_combination_count() << "\n";
        general_info_file << "t_max: " << t_max << "\n";
        general_info_file << "timeout: " << timeout << "\n";
        general_info_file << "max_threads: " << std::thread::hardware_concurrency() << "\n";
        general_info_file << "mechanism: " << par->model << "\n";
        general_info_file << "number_of_species: " << par->num_species << "\n";
        general_info_file << "number_of_reactions: " << par->num_reactions << "\n";
        general_info_file << "species: " << ::to_string((std::string*)par->species_names.data(), par->num_species) << "\n\n";
        general_info_file << parameter_combinator.to_string(true);
        general_info_file.close();
    }

    // save JSON file
    std::filesystem::path json_file_path = save_folder_path / "bruteforce_parameter_study_settings.json";
    std::ofstream json_file(json_file_path);
    if (!json_file.is_open())
    {
        LOG_ERROR("Failed to open file: " + json_file_path.string());
    } else {
        ordered_json j = parameter_combinator.to_json();
        j["info"] = ordered_json::object({
            {"version", VERSION},
            {"datetime", Timer::current_time()},
            {"total_combination_count", parameter_combinator.get_total_combination_count()},
            {"t_max", t_max},
            {"timeout", timeout},
            {"max_threads", std::thread::hardware_concurrency()},
            {"mechanism", Parameters::mechanism_names.at(parameter_combinator.mechanism)},
            {"number_of_species", parameter_combinator.get_mechanism_parameters()->num_species},
            {"number_of_reactions", parameter_combinator.get_mechanism_parameters()->num_reactions},
            {"species", ::to_string((std::string*)parameter_combinator.get_mechanism_parameters()->species_names.data(), parameter_combinator.get_mechanism_parameters()->num_species)}
        });
        json_file << std::setw(4) << j << std::endl;
        json_file.close();
    }

    // open log file
    std::filesystem::path log_file_path = save_folder_path / "output.log";
    std::filesystem::path error_file_path = save_folder_path / "errors.log";
    ErrorHandler::set_log_file(error_file_path.string());
    this->output_log_file.open(log_file_path);
    if (!this->output_log_file.is_open())
    {
        LOG_ERROR("Failed to open file: " + log_file_path.string());
    }
}


ParameterStudy::~ParameterStudy()
{
    if (this->output_log_file.is_open())
    {
        this->output_log_file.close();
    }
}


void ParameterStudy::parameter_study_task(const bool print_output, const size_t thread_id)
{
    // setup ODE solver and ODE function
    OdeFun ode;
    const Parameters* par = this->parameter_combinator.get_mechanism_parameters();
    if(par == nullptr) return;
    OdeSolverCVODE solver(4 + par->num_species);

    // open csv file
    std::filesystem::path save_folder_path(save_folder);
    std::filesystem::path csv_file_path = save_folder_path / ("output_" + std::to_string(thread_id) + ".csv");
    std::ofstream csv_file(csv_file_path);
    if (!csv_file.is_open())
    {
        LOG_ERROR("Failed to open file: " + csv_file_path.string());
        return;
    }
    csv_file << SimulationData::csv_header << "\n";

    // Setup SUNDIALS Logger directory
    std::filesystem::path sundials_log_path = save_folder_path / "sundials_logs";
    std::filesystem::create_directories(sundials_log_path);
    std::filesystem::path cvode_log_file_path = sundials_log_path / (std::to_string(thread_id) + ".log");
    solver.set_log_file(cvode_log_file_path.string());

    // run parameter study
    while(true)
    {
        // simulation setup
        auto [success, cpar] = parameter_combinator.get_next_combination();
        if (!success) break;
        ode.init(cpar);

        // run simulation and postprocess
        OdeSolution sol = solver.solve(
            this->t_max,                // t_max [s]
            &ode,                       // ode_ptr
            this->timeout,              // timeout [s]
            false                       // save solution
        );
        SimulationData data(cpar, sol);

        // save and print data
        csv_file << data.to_csv() << "\n";
        if (data.sol.success())
        {
            std::unique_lock<std::mutex> lock(this->output_mutex);
            this->best_energy_demand = std::min(this->best_energy_demand, data.energy_demand);
            lock.unlock();
        }

        if (data.sol.success())
        {
            this->successful_simulations.fetch_add(1, std::memory_order_relaxed); ;
        }
        this->total_simulations.fetch_add(1, std::memory_order_relaxed); ;

        if (this->output_log_file.is_open())
            this->output_log_file << data.to_small_string(parameter_combinator, best_energy_demand, false) << "\n";

        if (print_output)
        {
            std::unique_lock<std::mutex> lock(this->output_mutex);
            std::cout << data.to_small_string(parameter_combinator, best_energy_demand, true) << "\n";
            lock.unlock();
        }
    }

    csv_file.close();
}


void ParameterStudy::run(const size_t num_threads, const bool print_output)
{
    // Checks
    if (this->save_folder.empty()) return;
    if (num_threads == 0)
    {
        LOG_ERROR("Number of threads must be greater than zero.");
        return;
    }
    if (num_threads > std::thread::hardware_concurrency())
    {
        LOG_ERROR("Number of threads (" + std::to_string(num_threads) + ") must be less than or equal to the number of hardware threads (" + std::to_string(std::thread::hardware_concurrency()) + ")");
        return;
    }
    ErrorHandler::print_when_log = print_output;
    std::cout << "Running parameter study with " << num_threads << " threads..." << std::endl;

    // Run tasks
    Timer timer;
    timer.start();
    std::vector<std::thread> threads(num_threads);
    for (size_t i = 0; i < num_threads; i++)
    {
        threads[i] = std::thread(&ParameterStudy::parameter_study_task, this, print_output, i);
    }

    for (size_t i = 0; i < num_threads; i++)
    {
        threads[i].join();
    }

    // Print summary
    double runtime = timer.lap();
    std::cout << "\n\n";
    std::stringstream ss;
    ss << "Successful simulations: " << this->successful_simulations << "/" << this->total_simulations << " (" << std::setprecision(4) << 100.0 * this->successful_simulations / this->total_simulations << " %)" << std::endl;
    ss << "Total runtime: " << Timer::format_time(runtime) << std::endl;
    ss << "Average runtime per combination: " << Timer::format_time(runtime / parameter_combinator.get_total_combination_count()) << std::endl;
    std::cout << ss.str();

    // Save summary in txt
    std::filesystem::path save_folder_path(this->save_folder);
    std::filesystem::path summary_file_path = save_folder_path / "summary.txt";
    std::ofstream summary_file(summary_file_path);
    if (!summary_file.is_open())
    {
        LOG_ERROR("Failed to open file: " + summary_file_path.string());
    } else {
        summary_file << ss.str();
        summary_file.close();
    }
}

#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>

#include "nlohmann/json.hpp"
#include "parameter_combinator.h"
#include "ode_fun.h"

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


LogRange::LogRange(double start, double end, size_t num_steps):
    Range(start, end, std::max(num_steps, (size_t)1))
{
    if (start * end == 0.0) {
        LOG_ERROR("LogRange start and end must not be 0.");
        this->start = this->end = 1.0;
    }
    if (start * end < 0.0) {
        LOG_ERROR("LogRange start and end must both be positive or both negative.");
        this->start = this->end = 1.0;
    }

    if (start < 0 && end < 0) {
        this->sign = -1;
        this->start = std::abs(start);
        this->end = std::abs(end);
    } else {
        this->sign = 1;
    }

    this->a = std::log(this->start);
    if (this->num_steps > 1)
        this->b = (std::log(this->end) - std::log(this->start)) / (this->num_steps - 1);
    else
        this->b = 0.0;
}


double LogRange::operator[](size_t i) const
{
    if (this->num_steps <= 1)
    {
        return this->start;
    }
    if (i >= this->num_steps)
    {
        return this->end;
    }

    return this->sign * std::exp(this->a + this->b * i);
}


std::string LogRange::to_string() const
{
    return "LogRange(" + std::to_string(this->start) + ", " + std::to_string(this->end) + ", " + std::to_string(this->num_steps) + ")";
}


ordered_json LogRange::to_json() const
{
    return {
        {"type", "LogRange"},
        {"start", this->sign * this->start},
        {"end", this->sign * this->end},
        {"num_steps", this->num_steps}
    };
}


GeomRange::GeomRange(double start, double end, size_t num_steps, double q):
    Range(start, end, std::max(num_steps, (size_t)1)),
    q(q)
{
    if (this->num_steps == 1)
    {
        this->a_1 = end - start;
    }
    else
    {
        this->a_1 = (end - start) * (1.0 - q) / (1.0 - std::pow(q, this->num_steps - 1));
    }
}


double GeomRange::operator[](size_t i) const
{
    if (this->num_steps <= 1)
    {
        return this->start;
    }
    if (i >= this->num_steps)
    {
        return this->end;
    }

    return this->start + this->a_1 * (1.0 - std::pow(this->q, i)) / (1.0 - this->q);
}


std::string GeomRange::to_string() const
{
    return "GeomRange(" + std::to_string(this->start) + ", " + std::to_string(this->end) + ", " + std::to_string(this->num_steps) + ", " + std::to_string(this->q) + ")";
}


ordered_json GeomRange::to_json() const
{
    return {
        {"type", "GeomRange"},
        {"start", this->start},
        {"end", this->end},
        {"num_steps", this->num_steps},
        {"q", this->q}
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
    else if (std::holds_alternative<LogRange>(range))
    {
        return std::make_unique<LogRange>(std::get<LogRange>(range));
    }
    else if (std::holds_alternative<GeomRange>(range))
    {
        return std::make_unique<GeomRange>(std::get<GeomRange>(range));
    }
    return nullptr;
}


void ParameterCombinator::init(const ParameterCombinator::Builder &builder)
{
    this->combination_ID =              0;
    this->total_combination_count =     1;
    this->mechanism =                   builder.mechanism;
    this->R_E =                         get_unique_ptr(builder.R_E);
    this->ratio =                       get_unique_ptr(builder.ratio);
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
    this->total_combination_count *= this->ratio->get_num_steps();
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
            ". Available options: Const(value), LinearRange(start, end, num_steps), LogRange(start, end, num_steps), GeomRange(start, end, num_steps, q). " + \
            "E.g.: {\"type\": \"LinearRange\", \"start\": 0.0, \"end\": 1.0, \"num_steps\": 10}"
        );
        return default_range;
    }
    if (!j.contains("type") || !j.at("type").is_string())
    {
        LOG_ERROR(
            "Invalid \"type\" in JSON object " + j.dump() + \
            ". Available options: Const(value), LinearRange(start, end, num_steps), LogRange(start, end, num_steps), GeomRange(start, end, num_steps, q). " + \
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
    else if (type == "logrange")
    {
        if (
            !j.contains("start") || !j.contains("end") || !j.contains("num_steps")
            || !j.at("start").is_number() || !j.at("end").is_number() || !j.at("num_steps").is_number()
        )
        {
            LOG_ERROR(
                "Invalid arguments for JSON object " + j.dump() + \
                ". For LogRange, you must provide start, end, and optionally base of type double. As well as num_steps of type integer. " + \
                "e.g.: {\"type\": \"LogRange\", \"start\": 0.0, \"end\": 1.0, \"num_steps\": 10}"
            );
            return default_range;
        }
        return LogRange(
            j.at("start").get<double>(),
            j.at("end").get<double>(),
            j.at("num_steps").get<size_t>()
        );
    }
    else if (type == "geomrange")
    {
        if (
            !j.contains("start") || !j.contains("end") || !j.contains("num_steps") || !j.contains("q")
            || !j.at("start").is_number() || !j.at("end").is_number() || !j.at("num_steps").is_number() || !j.at("q").is_number()
        )
        {
            LOG_ERROR(
                "Invalid arguments for JSON object " + j.dump() + \
                ". For GeomRange, you must provide start, end, and q of type double. As well as num_steps of type integer. " + \
                "e.g.: {\"type\": \"GeomRange\", \"start\": 0.0, \"end\": 1.0, \"num_steps\": 10, \"q\": 2.0}"
            );
            return default_range;
        }
        return GeomRange(
            j.at("start").get<double>(),
            j.at("end").get<double>(),
            j.at("num_steps").get<size_t>(),
            j.at("q").get<double>()
        );
    }
    else
    {
        LOG_ERROR("Invalid Range type: " + type + \
            ". Available options: Const(value), LinearRange(start, end, num_steps), LogRange(start, end, num_steps), GeomRange(start, end, num_steps, q). " + \
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
        builder.R_E =                       get_range                           (j, "R_E",                      builder.R_E);
        builder.ratio =                     get_range                           (j, "ratio",                    builder.ratio);
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
        else
        {
            LOG_ERROR(
                Error::severity::warning,
                Error::type::preprocess,
                std::string("Warning, key excitation_params is missing. Excitation params are: ") + Parameters::excitation_arg_names.at(builder.excitation_type)
            );
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
    ss << format_field_name << ".ratio" << " = " << format_range(this->ratio.get()) << "\n";
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
    j["ratio"] = this->ratio->to_json();
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
        .ratio = get_value(this->ratio),
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

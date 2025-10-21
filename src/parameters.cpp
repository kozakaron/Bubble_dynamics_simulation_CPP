#include <algorithm>
#include <sstream>
#include <fstream>
#include <cmath>

#include "nlohmann/json.hpp"
#include "common.h"
#include "parameters.h"

using ordered_json = nlohmann::ordered_json;


/*#define COPY_ARRAY(type, name, size) { \
    type* temp = size ? new type[size] : nullptr; \
    this->name = (const type*)temp; \
    std::copy((type*)(T::name), (type*)(T::name) + size, temp); \
}*/


void compute_interval_values_and_derivatives(
    const double T_interval, // T_low or T_high
    const double* a,         // NASA coefficients: {a1, a2, a3, a4, a5, a6, a7}
    double* values,          // values: {C_p_interval, H_interval, S_interval}
    double* derivatives      // derivatives: {dC_p_interval/dT, dH_interval/dT, dS_interval/dT}
)
{
    // Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
    // H/RT = a1 + a2*T /2 + a3*T^2 /3 + a4*T^3 /4 + a5*T^4 /5 + a6/T
    // S/R  = a1*lnT + a2*T + a3*T^2 /2 + a4*T^3 /3 + a5*T^4 /4 + a7
    // dCp/dT = R * (a2 + 2*a3*T + 3*a4*T^2 + 4*a5*T^3)
    // dH/dT = R * (a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4)
    // dS/dT = R * (a1/T + a2 + a3*T + a4*T^2 + a5*T^3)

    const double T1 = T_interval;
    const double T2 = T1 * T1;
    const double T3 = T2 * T1;
    const double T4 = T2 * T2;
    const double T5 = T3 * T2;

    values[0] = Parameters::R_g * (a[0] + a[1]*T1 + a[2]*T2 + a[3]*T3 + a[4]*T4);   // C_p
    values[1] = Parameters::R_g * (a[0]*T1 + a[1]*T2/2.0 + a[2]*T3/3.0 + a[3]*T4/4.0 + a[4]*T5/5.0 + a[5]);   // H
    values[2] = Parameters::R_g * (a[0]*std::log(T1) + a[1]*T1 + a[2]*T2/2.0 + a[3]*T3/3.0 + a[4]*T4/4.0 + a[6]);   // S

    derivatives[0] = Parameters::R_g * (a[1] + 2.0*a[2]*T1 + 3.0*a[3]*T2 + 4.0*a[4]*T3);   // dC_p/dT
    derivatives[1] = Parameters::R_g * (a[0] + a[1]*T1 + a[2]*T2 + a[3]*T3 + a[4]*T4);   // dH/dT
    derivatives[2] = Parameters::R_g * (a[0]/T1 + a[1] + a[2]*T1 + a[3]*T2 + a[4]*T3);   // dS/dT
}


template<typename T>
T get_value(const nlohmann::ordered_json& j, const std::string& main_key, const std::string& sub_key="")
{
    constexpr bool is_floating_point = std::is_same_v<T, double> || std::is_same_v<T, float>;
    constexpr bool is_integer = std::is_same_v<T, int8_t> || std::is_same_v<T, int16_t> || std::is_same_v<T, int32_t> || std::is_same_v<T, int64_t> ||
                                std::is_same_v<T, uint8_t> || std::is_same_v<T, uint16_t> || std::is_same_v<T, uint32_t> || std::is_same_v<T, uint64_t>;
    constexpr bool is_string = std::is_same_v<T, std::string>;
    constexpr bool is_bool = std::is_same_v<T, bool>;
    constexpr bool is_float_vector = std::is_same_v<T, std::vector<double>> || std::is_same_v<T, std::vector<float>>;
    constexpr bool is_string_vector = std::is_same_v<T, std::vector<std::string>>;

    // check existence
    if (!j.contains(main_key))
    {
        LOG_ERROR(
            Error::severity::warning, Error::type::preprocess,
            "Key \"" + main_key + "\" not found in JSON object. Using default value. ");
        return T{};
    }

    if (!sub_key.empty() && !j.at(main_key).contains(sub_key))
    {
        LOG_ERROR(
            Error::severity::warning, Error::type::preprocess,
            "Key \"" + sub_key + "\" not found in JSON object. Using default value. ");
        return T{};
    }
    const nlohmann::ordered_json& target = sub_key.empty() ? j.at(main_key) : j.at(main_key).at(sub_key);

    // check type
    if (is_floating_point && !target.is_number_float())
    {
        LOG_ERROR(
            Error::severity::warning, Error::type::preprocess,
            "Expected floating point number for key \"" + sub_key + "\", instead found " + target.dump() + ". Using default value."
        );
        return T{};
    }
    else if (is_integer && !target.is_number_integer())
    {
        LOG_ERROR(
            Error::severity::warning, Error::type::preprocess,
            "Expected integer number for key \"" + sub_key + "\", instead found " + target.dump() + ". Using default value."
        );
        return T{};
    }
    else if (is_string && !target.is_string())
    {
        LOG_ERROR(
            Error::severity::warning, Error::type::preprocess,
            "Expected string for key \"" + sub_key + "\", instead found " + target.dump() + ". Using default value."
        );
        return T{};
    }
    else if (is_bool && !target.is_boolean())
    {
        LOG_ERROR(
            Error::severity::warning, Error::type::preprocess,
            "Expected boolean for key \"" + sub_key + "\", instead found " + target.dump() + ". Using default value."
        );
        return T{};
    }
    else if ((is_float_vector || is_string_vector) && !target.is_array())
    {
        LOG_ERROR(
            Error::severity::warning, Error::type::preprocess,
            "Expected JSON array for key \"" + sub_key + "\", instead found " + target.dump() + ". Using default value."
        );
        return T{};
    }

    // check element types for vectors
    if (is_float_vector || is_string_vector)
    {
        for (const auto& element : target)
        {
            if (is_float_vector && !element.is_number_float())
            {
                LOG_ERROR(
                    Error::severity::warning, Error::type::preprocess,
                    "Expected floating point number in array for key \"" + sub_key + "\", instead found " + element.dump() + ". Using default value."
                );
                return T{};
            }
            else if (is_string_vector && !element.is_string())
            {
                LOG_ERROR(
                    Error::severity::warning, Error::type::preprocess,
                    "Expected string in array for key \"" + sub_key + "\", instead found " + element.dump() + ". Using default value."
                );
                return T{};
            }
        }
    }

    // get element
    return target.get<T>();
}


template<typename T>
void copy_array_1d(const nlohmann::ordered_json& j, const std::string& main_key, const std::string& sub_key, const T*& array, const index_t size)
{
    std::vector<T> input = get_value<std::vector<T>>(j, main_key, sub_key);
    if (size == 0 || input.size() != size)
    {
        array = nullptr;
        return;
    }
    T* temp_array = new T[size];
    for (index_t i = 0; i < size; i++)
    {
        temp_array[i] = input.at(i);
    }
    array = (const T*)temp_array;
}


template<typename T>
void copy_array_2d(const nlohmann::ordered_json& j, const std::string& main_key, const std::string& sub_key, const T*& array, const index_t rows, const index_t cols)
{
    std::vector<std::vector<T>> input = get_value<std::vector<std::vector<T>>>(j, main_key, sub_key);
    if (rows == 0 || cols == 0 || input.size() != rows || input[0].size() != cols)
    {
        array = nullptr;
        return;
    }
    T* temp_array = new T[rows*cols];
    for (index_t i = 0; i < rows; i++)
    {
        if (input.at(i).size() != cols)
        {
            delete[] temp_array;
            array = nullptr;
            return;
        }

        T* row_ptr = temp_array + i * cols;
        for (index_t j = 0; j < cols; j++)
        {
            row_ptr[j] = input.at(i).at(j);
        }
    }
    array = (const T*)temp_array;
}

Parameters::Parameters(const nlohmann::ordered_json& j):
    _species(),
    model(                          get_value<std::string>(j, "model")),
    num_elements(                   get_value<index_t>(j, "species", "num_elements")),
    num_species(                    get_value<index_t>(j, "species", "num_species")),
    index_of_water(                 get_value<index_t>(j, "species", "index_of_water")),
    invalid_index(                  get_value<index_t>(j, "species", "invalid_index")),
    num_reactions(                  get_value<index_t>(j, "arrhenius_parameters", "num_reactions")),
    num_max_species_per_reaction(   get_value<index_t>(j, "arrhenius_parameters", "num_max_species_per_reaction")),
    num_third_body_reactions(       get_value<index_t>(j, "third_body_reactions", "num_third_body_reactions")),
    num_irreversible_reactions(     get_value<index_t>(j, "irreversible_reactions", "num_irreversible_reactions")),
    num_falloff_reactions(          get_value<index_t>(j, "falloff_reactions", "num_falloff_reactions")),
    num_lindemann_reactions(        get_value<index_t>(j, "falloff_reactions", "num_lindemann_reactions")),
    num_troe_reactions(             get_value<index_t>(j, "falloff_reactions", "num_troe_reactions")),
    num_sri_reactions(              get_value<index_t>(j, "falloff_reactions", "num_sri_reactions")),
    num_plog_reactions(             get_value<index_t>(j, "plog_reactions", "num_plog_reactions")),
    num_plog_levels(                get_value<index_t>(j, "plog_reactions", "num_plog_levels")),
    molar_weights(nullptr),
    thermal_conductivities(nullptr),
    temp_ranges(nullptr),
    a_low(nullptr),
    a_high(nullptr),
    interval_values(nullptr),
    interval_derivatives(nullptr),
    b(nullptr),
    ln_A(nullptr),
    E_over_R(nullptr),
    reaction_order(nullptr),
    nu_indexes(nullptr),
    nu_forward(nullptr),
    nu_backward(nullptr),
    nu(nullptr),
    sum_nu(nullptr),
    third_body_reaction_indexes(nullptr),
    is_falloff_reaction(nullptr),
    third_body_efficiencies(nullptr),
    irreversible_reaction_indexes(nullptr),
    falloff_reaction_indexes(nullptr),
    pressure_dependent_reac_types(nullptr),
    is_third_body_indexes(nullptr),
    falloff_parameters(nullptr),
    troe_parameters(nullptr),
    sri_parameters(nullptr),
    plog_reaction_indexes(nullptr),
    plog_seperators(nullptr),
    plog_parameters(nullptr)
{
    std::vector<std::string> main_keys = {
        "model",
        "species",
        "thermodynamics",
        "arrhenius_parameters",
        "third_body_reactions",
        "irreversible_reactions",
        "falloff_reactions",
        "plog_reactions"
    };
    for (const auto& key: main_keys)
    {
        if (!j.contains(key))
        {
            LOG_ERROR(
                Error::severity::error, Error::type::preprocess,
                "JSON file does not contain '" + key + "' key. It might be the wrong json: " + j.dump()
            );
            return;
        }
    }

    // Fill simple arrays
    copy_array_1d<double> (j, "species", "molar_weights", molar_weights, num_species);
    copy_array_1d<double> (j, "species", "thermal_conductivities", thermal_conductivities, num_species);
    copy_array_2d<double> (j, "thermodynamics", "temp_ranges", temp_ranges, num_species, 3);
    copy_array_2d<double> (j, "thermodynamics", "a_low", a_low, num_species, NASA_order + 2);
    copy_array_2d<double> (j, "thermodynamics", "a_high", a_high, num_species, NASA_order + 2);
    copy_array_1d<index_t>(j, "third_body_reactions", "third_body_reaction_indexes", third_body_reaction_indexes, num_third_body_reactions);
    copy_array_1d<bool>   (j, "third_body_reactions", "is_falloff_reaction", is_falloff_reaction, num_third_body_reactions);
    copy_array_2d<double> (j, "third_body_reactions", "third_body_efficiencies", third_body_efficiencies, num_third_body_reactions, num_species);
    copy_array_1d<index_t>(j, "irreversible_reactions", "irreversible_reaction_indexes", irreversible_reaction_indexes, num_irreversible_reactions);
    copy_array_1d<index_t>(j, "falloff_reactions", "falloff_reaction_indexes", falloff_reaction_indexes, num_falloff_reactions);
    copy_array_1d<index_t>(j, "falloff_reactions", "is_third_body_indexes", is_third_body_indexes, num_falloff_reactions);
    copy_array_2d<double> (j, "falloff_reactions", "falloff_parameters", falloff_parameters, num_falloff_reactions, 3);
    copy_array_2d<double> (j, "falloff_reactions", "troe_parameters", troe_parameters, num_troe_reactions, 4);
    copy_array_2d<double> (j, "falloff_reactions", "sri_parameters", sri_parameters, num_sri_reactions, 5);
    copy_array_1d<index_t>(j, "plog_reactions", "plog_reaction_indexes", plog_reaction_indexes, num_plog_reactions);
    copy_array_1d<index_t>(j, "plog_reactions", "plog_seperators", plog_seperators, num_plog_levels + 1);
    copy_array_2d<double> (j, "plog_reactions", "plog_parameters", plog_parameters, num_plog_levels, 4);
    
    // Species map and vector
    species_names = get_value<std::vector<std::string>>(j, "species", "species_names");
    for (index_t i = 0; i < species_names.size(); i++)
    {
        _species[species_names.at(i)] = i;
    }

    // Arrhenius parameters and reaction order
    // TODO

    // Stoichiometric coefficients
    // TODO

    // Falloff reaction types
    // TODO

}

// Removes "//" comments from a line
// Does not support the placement of "//" in strings
std::string strip_comments(const std::string& line) {
    size_t pos = line.find("//");
    return (pos != std::string::npos) ? line.substr(0, pos) : line;
}


ordered_json get_json(const std::string& json_path)
{
    // Open JSON file
    std::ifstream input_file(json_path);
    if (!input_file.is_open())
    {
        LOG_ERROR(
            Error::severity::error,
            Error::type::preprocess,
            "Could not open JSON file: " + json_path
        );
        return {};
    }

    // Read JSON data
    std::ostringstream cleaned;
    std::string line;
    ordered_json j;
    while (std::getline(input_file, line)) {
        cleaned << strip_comments(line) << '\n';
    }
    try
    {
        j = ordered_json::parse(cleaned.str());
    }
    catch (const std::exception& e)
    {
        LOG_ERROR(
            Error::severity::error,
            Error::type::preprocess,
            "Error parsing JSON file: " + std::string(e.what())
        );
        return {};
    }
    input_file.close();

    // Parse JSON data
    if (!j.contains("model"))
    {
        LOG_ERROR(
            Error::severity::error,
            Error::type::preprocess,
            "JSON file does not contain 'model' key. It might be the wrong json: " + j.dump()
        );
        return {};
    }
    
    return j;
}


Parameters::Parameters(const std::string& json_path):
    Parameters(get_json(json_path))
{ }


/*
template <typename T>
Parameters::Parameters(T dummy):
    _elements(),
    _species(),
    model(T::model),
    input_file(T::input_file),
    num_elements(T::num_elements),
    num_species(T::num_species),
    index_of_water(T::index_of_water),
    invalid_index(T::invalid_index),
    species_names(),
    NASA_order(T::NASA_order),
    num_reactions(T::num_reactions),
    num_max_species_per_reaction(T::num_max_species_per_reaction),
    num_third_body_reactions(T::num_third_body_reactions),
    num_irreversible_reactions(T::num_irreversible_reactions),
    num_falloff_reactions(T::num_falloff_reactions),
    num_lindemann_reactions(T::num_lindemann_reactions),
    num_troe_reactions(T::num_troe_reactions),
    num_sri_reactions(T::num_sri_reactions),
    num_plog_reactions(T::num_plog_reactions),
    num_plog_levels(T::num_plog_levels)
{
    (void)dummy;
    
    // Physical constants
    COPY_ARRAY(double, molar_weights, T::num_species);
    COPY_ARRAY(double, thermal_conductivities, T::num_species);

    // NASA polynomials
    COPY_ARRAY(double, temp_ranges, T::num_species*3);
    COPY_ARRAY(double, a_low, T::num_species*(T::NASA_order+2));
    COPY_ARRAY(double, a_high, T::num_species*(T::NASA_order+2));
    double* interval_values_temp = new double[T::num_species*3];
    double* interval_derivatives_temp = new double[T::num_species*3];
    this->interval_values = (const double*)interval_values_temp;
    this->interval_derivatives = (const double*)interval_derivatives_temp;
    for(index_t i = 0; i < T::num_species; i++)
    {
        compute_interval_values_and_derivatives(
            T::temp_ranges[i][1],    // T_high
            &(T::a_high[i][0]),
            &(interval_values_temp[i*3]),
            &(interval_derivatives_temp[i*3])
        );
    }

    // Arrhenius parameters
    COPY_ARRAY(double, b, T::num_reactions);
    COPY_ARRAY(double, ln_A, T::num_reactions);
    COPY_ARRAY(double, E_over_R, T::num_reactions);

    // Reaction order and stoichiometric coefficients
    COPY_ARRAY(index_t, reaction_order, T::num_reactions);
    COPY_ARRAY(index_t, nu_indexes, T::num_reactions*T::num_max_species_per_reaction);
    for(index_t i = 0; i < T::num_elements; i++)
    {
        _elements[T::elements[i].first] = T::elements[i].second;
        elements_names.push_back(T::elements[i].first);
    }
    for(index_t i = 0; i < T::num_species; i++)
    {
        _species[T::species[i].first] = T::species[i].second;
        species_names.push_back(T::species[i].first);
    }
    stoich_t *temp_nu_forward = new stoich_t[T::num_reactions*T::num_max_species_per_reaction];
    stoich_t *temp_nu_backward = new stoich_t[T::num_reactions*T::num_max_species_per_reaction];
    stoich_t *temp_nu = new stoich_t[T::num_reactions*T::num_max_species_per_reaction];
    stoich_t *temp_sum_nu = new stoich_t[T::num_reactions];
    for(index_t i = 0; i < T::num_reactions; i++)
    {
        std::copy(
            (stoich_t*)T::nu[i][0],
            (stoich_t*)T::nu[i][0] + T::num_max_species_per_reaction,
            temp_nu_forward + i*T::num_max_species_per_reaction
        );
        std::copy(
            (stoich_t*)T::nu[i][1],
            (stoich_t*)T::nu[i][1] + T::num_max_species_per_reaction,
            temp_nu_backward + i*T::num_max_species_per_reaction
        );
        std::copy(
            (stoich_t*)T::nu[i][2],
            (stoich_t*)T::nu[i][2] + T::num_max_species_per_reaction,
            temp_nu + i*T::num_max_species_per_reaction
        );
        temp_sum_nu[i] = 0;
        for(index_t j = 0; j < T::num_max_species_per_reaction; j++)
        {
            temp_sum_nu[i] += T::nu[i][2][j];
        }
    }
    this->nu_forward = (const stoich_t*)temp_nu_forward;
    this->nu_backward = (const stoich_t*)temp_nu_backward;
    this->nu = (const stoich_t*)temp_nu;
    this->sum_nu = (const stoich_t*)temp_sum_nu;

    // Third body and fallof  pressure-dependent reactions
    COPY_ARRAY(index_t, third_body_reaction_indexes, T::num_third_body_reactions);
    COPY_ARRAY(bool, is_falloff_reaction, T::num_third_body_reactions);
    COPY_ARRAY(double, third_body_efficiencies, T::num_third_body_reactions*T::num_species);
    COPY_ARRAY(index_t, irreversible_reaction_indexes, T::num_irreversible_reactions);
    COPY_ARRAY(index_t, falloff_reaction_indexes, T::num_falloff_reactions);
    COPY_ARRAY(Parameters::reac_type, pressure_dependent_reac_types, T::num_falloff_reactions);
    COPY_ARRAY(index_t, is_third_body_indexes, T::num_falloff_reactions);
    COPY_ARRAY(double, falloff_parameters, T::num_falloff_reactions*3);
    COPY_ARRAY(double, troe_parameters, T::num_troe_reactions*4);
    COPY_ARRAY(double, sri_parameters, T::num_sri_reactions*5);

    // PLOG reactions
    COPY_ARRAY(index_t, plog_reaction_indexes, T::num_plog_reactions);
    if (T::num_plog_reactions > 0)
    {
        COPY_ARRAY(index_t, plog_seperators, T::num_plog_reactions+1);
    } else {
        COPY_ARRAY(index_t, plog_seperators, 2);
    }
   COPY_ARRAY(double, plog_parameters, T::num_plog_levels*4);
}*/


Parameters::~Parameters()
{
    if (molar_weights != nullptr) delete[] molar_weights;
    if (thermal_conductivities != nullptr) delete[] thermal_conductivities;
    if (temp_ranges != nullptr) delete[] temp_ranges;
    if (a_low != nullptr) delete[] a_low;
    if (a_high != nullptr) delete[] a_high;
    if (interval_values != nullptr) delete[] interval_values;
    if (interval_derivatives != nullptr) delete[] interval_derivatives;
    if (b != nullptr) delete[] b;
    if (ln_A != nullptr) delete[] ln_A;
    if (E_over_R != nullptr) delete[] E_over_R;
    if (nu_indexes != nullptr) delete[] nu_indexes;
    if (nu_forward != nullptr) delete[] nu_forward;
    if (nu_backward != nullptr) delete[] nu_backward;
    if (reaction_order != nullptr) delete[] reaction_order;
    if (nu != nullptr) delete[] nu;
    if (sum_nu != nullptr) delete[] sum_nu;
    if (third_body_reaction_indexes != nullptr) delete[] third_body_reaction_indexes;
    if (is_falloff_reaction != nullptr) delete[] is_falloff_reaction;
    if (third_body_efficiencies != nullptr) delete[] third_body_efficiencies;
    if (irreversible_reaction_indexes != nullptr) delete[] irreversible_reaction_indexes;
    if (falloff_reaction_indexes != nullptr) delete[] falloff_reaction_indexes;
    if (pressure_dependent_reac_types != nullptr) delete[] pressure_dependent_reac_types;
    if (is_third_body_indexes != nullptr) delete[] is_third_body_indexes;
    if (falloff_parameters != nullptr) delete[] falloff_parameters;
    if (troe_parameters != nullptr) delete[] troe_parameters;
    if (sri_parameters != nullptr) delete[] sri_parameters;
    if (plog_reaction_indexes != nullptr) delete[] plog_reaction_indexes;
    if (plog_seperators != nullptr) delete[] plog_seperators;
    if (plog_parameters != nullptr) delete[] plog_parameters;
}

const Parameters *Parameters::get_parameters(const std::string& mech)
{
    /*switch (mech)
    {
        case mechanism::chemkin_ar_he:
            return &Parameters::chemkin_ar_he_params;
        case mechanism::chemkin_kaust2023_n2:
            return &Parameters::chemkin_kaust2023_n2_params;
        case mechanism::chemkin_kaust2023_n2_without_o:
            return &Parameters::chemkin_kaust2023_n2_without_o_params;
        case mechanism::chemkin_otomo2018_without_o:
            return &Parameters::chemkin_otomo2018_without_o_params;
        case mechanism::chemkin_otomo2018:
            return &Parameters::chemkin_otomo2018_params;
        case mechanism::chemkin_elte2016_ethanol:
            return &Parameters::chemkin_elte2016_ethanol_params;
        case mechanism::chemkin_elte2016_syngas:
            return &Parameters::chemkin_elte2016_syngas_params;
        case mechanism::chemkin_elte2017_methanol:
            return &Parameters::chemkin_elte2017_methanol_params;
        default:
            LOG_ERROR("Unknown mechanisms: " + std::to_string(mech));
            return nullptr;
    }*/ return nullptr;
}


index_t Parameters::get_species(std::string name) const
{
    std::transform(name.begin(), name.end(), name.begin(), ::toupper);
    auto it = _species.find(name);
    if (it == _species.end())
    {
        std::stringstream ss;
        ss << "Species \"" << name << "\" not found in " << this->model << ". Valid species are: ";
        ss << ::to_string((std::string*)Parameters::species_names.data(), Parameters::species_names.size());
        LOG_ERROR(Error::severity::error, Error::type::preprocess, ss.str());
        return invalid_index;
    }
    return it->second;
}


/*Parameters::mechanism Parameters::string_to_mechanism(std::string& mechanism_str) {
    std::transform(mechanism_str.begin(), mechanism_str.end(), mechanism_str.begin(), ::tolower);
    std::replace(mechanism_str.begin(), mechanism_str.end(), ' ', '_');
    auto it = std::find(Parameters::mechanism_names.begin(), Parameters::mechanism_names.end(), mechanism_str);
    if (it != Parameters::mechanism_names.end()) {
        return static_cast<Parameters::mechanism>(std::distance(Parameters::mechanism_names.begin(), it));
    }

    // Error handling
    std::stringstream ss;
    ss << "Invalid mechanism: " << mechanism_str << ". Valid options are: ";
    ss << ::to_string((char**)Parameters::mechanism_names.data(), Parameters::mechanism_names.size());
    LOG_ERROR(Error::severity::error, Error::type::preprocess, ss.str());
    return Parameters::mechanism::chemkin_ar_he; // Default to a known mechanism
}*/


Parameters::excitation Parameters::string_to_excitation(std::string excitation_str) {
    std::transform(excitation_str.begin(), excitation_str.end(), excitation_str.begin(), ::tolower);
    std::replace(excitation_str.begin(), excitation_str.end(), ' ', '_');
    auto it = std::find(Parameters::excitation_names.begin(), Parameters::excitation_names.end(), excitation_str);
    if (it != Parameters::excitation_names.end()) {
        return static_cast<Parameters::excitation>(std::distance(Parameters::excitation_names.begin(), it));
    }

    // Error handling
    std::stringstream ss;
    ss << "Invalid excitation type: " << excitation_str << ". Valid options are: ";
    ss << ::to_string((char**)Parameters::excitation_names.data(), Parameters::excitation_names.size());
    LOG_ERROR(Error::severity::error, Error::type::preprocess, ss.str());
    return Parameters::excitation::sin_impulse; // Default to a known excitation type
}
#include <algorithm>
#include <sstream>
#include <fstream>
#include <cmath>
#include <filesystem>

#include "nlohmann/json.hpp"
#include "common.h"
#include "parameters.h"

using json = nlohmann::json;


std::vector<std::unique_ptr<Parameters>>  Parameters::_mechanisms{};


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


// Helper function to check if key and subkey exist in JSON
bool check_key_exists(
    const nlohmann::json& j,
    const std::string& main_key,
    const std::string& sub_key=""
)
{
    if (!j.contains(main_key))
    {
        LOG_ERROR(
            Error::severity::warning, Error::type::preprocess,
            "Key \"" + main_key + "\" not found in JSON object. ");
        return false;
    }

    if (!sub_key.empty() && !j.at(main_key).contains(sub_key))
    {
        LOG_ERROR(
            Error::severity::warning, Error::type::preprocess,
            "Subkey \"" + sub_key + "\" under key \"" + main_key + "\" not found in JSON object. ");
        return false;
    }

    return true;
}


// Helper function to check if key type matches expected type (for scalars and 1D vectors)
template<typename T>
bool check_key_type(const nlohmann::json& target)
{
    constexpr bool is_floating_point = std::is_same_v<T, double> || std::is_same_v<T, float>;
    constexpr bool is_integer = std::is_same_v<T, int8_t> || std::is_same_v<T, int16_t> || std::is_same_v<T, int32_t> || std::is_same_v<T, int64_t> ||
                                std::is_same_v<T, uint8_t> || std::is_same_v<T, uint16_t> || std::is_same_v<T, uint32_t> || std::is_same_v<T, uint64_t>;
    constexpr bool is_string = std::is_same_v<T, std::string>;
    constexpr bool is_bool = std::is_same_v<T, bool>;
    constexpr bool is_int_vector = std::is_same_v<T, std::vector<int8_t>> || std::is_same_v<T, std::vector<int16_t>> || std::is_same_v<T, std::vector<int32_t>> || std::is_same_v<T, std::vector<int64_t>> ||
                                   std::is_same_v<T, std::vector<uint8_t>> || std::is_same_v<T, std::vector<uint16_t>> || std::is_same_v<T, std::vector<uint32_t>> || std::is_same_v<T, std::vector<uint64_t>>;
    constexpr bool is_float_vector = std::is_same_v<T, std::vector<double>> || std::is_same_v<T, std::vector<float>>;
    constexpr bool is_string_vector = std::is_same_v<T, std::vector<std::string>>;
    constexpr bool is_bool_vector = std::is_same_v<T, std::vector<bool>>;
    constexpr bool is_vector = is_int_vector || is_float_vector || is_string_vector || is_bool_vector;

    std::string message = "";
    if (is_floating_point && !target.is_number_float())
    {
        message = "Expected floating point number for key, instead found " + target.dump() + ". Using default value.";
    }
    else if (is_integer && !target.is_number_integer())
    {
        message = "Expected integer number for key, instead found " + target.dump() + ". Using default value.";
    }
    else if (is_string && !target.is_string())
    {
        message = "Expected string for key, instead found " + target.dump() + ". Using default value.";
    }
    else if (is_bool && !target.is_boolean())
    {
        message = "Expected boolean for key, instead found " + target.dump() + ". Using default value.";
    }
    else if (is_vector && !target.is_array())
    {
        message = "Expected JSON array for key, instead found " + target.dump() + ". Using default value.";
    }

    if (!message.empty())
    {
        LOG_ERROR(Error::severity::warning, Error::type::preprocess, message);
        return false;
    }
    
    return true;
}


// Helper function to get value from JSON with type checking (for scalars and 1D vectors)
template<typename T>
T get_value(const nlohmann::json& j, const std::string& main_key, const std::string& sub_key="")
{
    constexpr bool is_int_vector = std::is_same_v<T, std::vector<int8_t>> || std::is_same_v<T, std::vector<int16_t>> || std::is_same_v<T, std::vector<int32_t>> || std::is_same_v<T, std::vector<int64_t>> ||
                                   std::is_same_v<T, std::vector<uint8_t>> || std::is_same_v<T, std::vector<uint16_t>> || std::is_same_v<T, std::vector<uint32_t>> || std::is_same_v<T, std::vector<uint64_t>>;
    constexpr bool is_float_vector = std::is_same_v<T, std::vector<double>> || std::is_same_v<T, std::vector<float>>;
    constexpr bool is_string_vector = std::is_same_v<T, std::vector<std::string>>;
    constexpr bool is_bool_vector = std::is_same_v<T, std::vector<bool>>;
    constexpr bool is_vector = is_int_vector || is_float_vector || is_string_vector || is_bool_vector;

    // check existence
    if (!check_key_exists(j, main_key, sub_key)) return T{};
    const nlohmann::json& target = sub_key.empty() ? j.at(main_key) : j.at(main_key).at(sub_key);

    // check type
    if (!check_key_type<T>(target)) return T{};

    // check element types for vectors
    if (is_vector)
    {
        for (const auto& element : target)
        {
            std::string message = "";
            if (is_int_vector && !element.is_number_integer())
            {
                message = "Expected integer number in array for key \"" + sub_key + "\", instead found " + element.dump() + ". Using default value.";
            }
            else if (is_float_vector && !element.is_number_float())
            {
                message = "Expected floating point number in array for key \"" + sub_key + "\", instead found " + element.dump() + ". Using default value.";
            }
            else if (is_string_vector && !element.is_string())
            {
                message = "Expected string in array for key \"" + sub_key + "\", instead found " + element.dump() + ". Using default value.";
            }
            else if (is_bool_vector && !element.is_boolean())
            {
                message = "Expected boolean in array for key \"" + sub_key + "\", instead found " + element.dump() + ". Using default value.";
            }

            if (!message.empty())
            {
                LOG_ERROR(Error::severity::warning, Error::type::preprocess, message);
                return T{};
            }
        }
    }

    // get element
    return target.get<T>();
}


// Helper function to copy 1D arrays from JSON to dynamically allocated C-style arrays
template<typename T>
void copy_array_1d(const nlohmann::json& j, const std::string& main_key, const std::string& sub_key, const T*& array, const index_t size)
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


// Helper function to copy 2D arrays from JSON to dynamically allocated C-style arrays
// NOTE: Element types in 2D arrays are not checked
template<typename T>
void copy_array_2d(const nlohmann::json& j, const std::string& main_key, const std::string& sub_key, const T*& array, const index_t rows, const index_t cols)
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


Parameters::Parameters(const nlohmann::json& j):
    _species(),
    mechanism_name(                 get_value<std::string>(j, "model")),
    num_elements(                   get_value<index_t>(j, "species", "num_elements")),
    num_species(                    get_value<index_t>(j, "species", "num_species")),
    index_of_water(                 get_value<index_t>(j, "species", "index_of_water")),
    invalid_index(                  get_value<index_t>(j, "species", "invalid_index")),
    W(                              nullptr),
    thermal_conductivities(         nullptr),
    sqrt_van_der_waals_a(           nullptr),
    van_der_waals_b(                nullptr),
    temp_ranges(                    nullptr),
    a_low(                          nullptr),
    a_high(                         nullptr),
    interval_values(                nullptr),
    interval_derivatives(           nullptr),
    num_reactions(                  get_value<index_t>(j, "arrhenius_parameters", "num_reactions")),
    arrhenius_parameters(           nullptr),
    reaction_order(                 nullptr),
    num_max_species_per_reaction(   get_value<index_t>(j, "arrhenius_parameters", "num_max_species_per_reaction")),
    nu_indexes(                     nullptr),
    nu_forward(                     nullptr),
    nu_backward(                    nullptr),
    nu(                             nullptr),
    sum_nu(                         nullptr),
    num_third_body_reactions(       get_value<index_t>(j, "third_body_reactions", "num_third_body_reactions")),
    third_body_reaction_indexes(    nullptr),
    is_falloff_reaction(            nullptr),
    alpha(                          nullptr),
    num_irreversible_reactions(     get_value<index_t>(j, "irreversible_reactions", "num_irreversible_reactions")),
    irreversible_reaction_indexes(  nullptr),
    num_falloff_reactions(          get_value<index_t>(j, "falloff_reactions", "num_falloff_reactions")),
    num_lindemann_reactions(        get_value<index_t>(j, "falloff_reactions", "num_lindemann_reactions")),
    num_troe_reactions(             get_value<index_t>(j, "falloff_reactions", "num_troe_reactions")),
    num_sri_reactions(              get_value<index_t>(j, "falloff_reactions", "num_sri_reactions")),
    falloff_reaction_indexes(       nullptr),
    falloff_reaction_types(         nullptr),
    is_third_body_indexes(          nullptr),
    falloff_parameters(             nullptr),
    troe_parameters(                nullptr),
    sri_parameters(                 nullptr),
    num_plog_reactions(             get_value<index_t>(j, "plog_reactions", "num_plog_reactions")),
    num_plog_levels(                get_value<index_t>(j, "plog_reactions", "num_plog_levels")),
    plog_reaction_indexes(          nullptr),
    plog_seperators(                nullptr),
    plog_parameters(                nullptr)
{
    // Fill simple arrays
    copy_array_1d<double> (j, "species", "molar_weights", W, num_species);
    copy_array_1d<double> (j, "species", "thermal_conductivities", thermal_conductivities, num_species);
    copy_array_1d<double> (j, "species", "van_der_waals_b", van_der_waals_b, num_species);
    copy_array_2d<double> (j, "thermodynamics", "temp_ranges", temp_ranges, num_species, 3);
    copy_array_2d<double> (j, "thermodynamics", "a_low", a_low, num_species, NASA_order + 2);
    copy_array_2d<double> (j, "thermodynamics", "a_high", a_high, num_species, NASA_order + 2);
    copy_array_1d<index_t>(j, "third_body_reactions", "third_body_reaction_indexes", third_body_reaction_indexes, num_third_body_reactions);
    copy_array_1d<bool>   (j, "third_body_reactions", "is_falloff_reaction", is_falloff_reaction, num_third_body_reactions);
    copy_array_2d<double> (j, "third_body_reactions", "third_body_efficiencies", alpha, num_third_body_reactions, num_species);
    copy_array_1d<index_t>(j, "irreversible_reactions", "irreversible_reaction_indexes", irreversible_reaction_indexes, num_irreversible_reactions);
    copy_array_1d<index_t>(j, "falloff_reactions", "falloff_reaction_indexes", falloff_reaction_indexes, num_falloff_reactions);
    copy_array_1d<index_t>(j, "falloff_reactions", "is_third_body_indexes", is_third_body_indexes, num_falloff_reactions);
    copy_array_2d<double> (j, "falloff_reactions", "falloff_parameters", falloff_parameters, num_falloff_reactions, 3);
    copy_array_2d<double> (j, "falloff_reactions", "troe_parameters", troe_parameters, num_troe_reactions, 4);
    copy_array_2d<double> (j, "falloff_reactions", "sri_parameters", sri_parameters, num_sri_reactions, 5);
    copy_array_1d<index_t>(j, "plog_reactions", "plog_reaction_indexes", plog_reaction_indexes, num_plog_reactions);
    copy_array_1d<index_t>(j, "plog_reactions", "plog_seperators", plog_seperators, num_plog_reactions + 1);
    copy_array_2d<double> (j, "plog_reactions", "plog_parameters", plog_parameters, num_plog_levels, 4);
    
    try
    {
    // Species map and vector
        species_names = get_value<std::vector<std::string>>(j, "species", "species_names");
        for (index_t i = 0; i < species_names.size(); i++)
        {
            _species[species_names.at(i)] = i;
        }

    // Sqare root of van der Waals a parameter
        if (!check_key_exists(j, "species", "van_der_waals_a")) return;
        const nlohmann::json& vdw_a_array = j.at("species").at("van_der_waals_a");
        double* sqrt_van_der_waals_a_temp = new double[num_species];
        for (index_t i = 0; i < num_species; i++)
        {
            double vdw_a = vdw_a_array.at(i).get<double>();
            sqrt_van_der_waals_a_temp[i] = std::sqrt(vdw_a);
        }
        sqrt_van_der_waals_a = (const double*)sqrt_van_der_waals_a_temp;


    // Compute NASA interval values and derivatives
        double* interval_values_temp      = new double[num_species * 3];
        double* interval_derivatives_temp = new double[num_species * 3];
        for(index_t i = 0; i < num_species; i++)
        {
            compute_interval_values_and_derivatives(
                temp_ranges[3*i + 2],    // T_high
                &(a_high[i*7]),
                &(interval_values_temp[i*3]),
                &(interval_derivatives_temp[i*3])
            );
        }
        interval_values      = (const double*)interval_values_temp;
        interval_derivatives = (const double*)interval_derivatives_temp;


    // Arrhenius parameters and reaction order
        if (!check_key_exists(j, "arrhenius_parameters", "arrhenius_params_with_orders")) return;
        const nlohmann::json& arrhenius = j.at("arrhenius_parameters").at("arrhenius_params_with_orders");
        if (arrhenius.size() != num_reactions) return;

        double* arrhenius_parameters_temp = new double[num_reactions * 3];
        index_t* reaction_order_temp      = new index_t[num_reactions];
        for (size_t i = 0; i < arrhenius.size(); i++)
        {
            const auto& reaction = arrhenius.at(i);
            arrhenius_parameters_temp[i * 3 + 0] = reaction.at(0).at(0).get<double>();
            arrhenius_parameters_temp[i * 3 + 1] = reaction.at(0).at(1).get<double>();
            arrhenius_parameters_temp[i * 3 + 2] = reaction.at(0).at(2).get<double>();
            reaction_order_temp[i] = reaction.at(1).get<index_t>();
        }
        arrhenius_parameters = (const double*)arrhenius_parameters_temp;
        reaction_order = (const index_t*)reaction_order_temp;


    // Stoichiometric coefficients
        if (!check_key_exists(j, "arrhenius_parameters", "stochiometric_coeffs")) return;
        const nlohmann::json& stoichiometric = j.at("arrhenius_parameters").at("stochiometric_coeffs");
        if (stoichiometric.size() != num_reactions) return;

        index_t* nu_indexes_temp   = new index_t[num_reactions * num_max_species_per_reaction];
        stoich_t* nu_forward_temp  = new stoich_t[num_reactions * num_max_species_per_reaction];
        stoich_t* nu_backward_temp = new stoich_t[num_reactions * num_max_species_per_reaction];
        stoich_t* nu_temp          = new stoich_t[num_reactions * num_max_species_per_reaction];
        stoich_t* sum_nu_temp      = new stoich_t[num_reactions];
        for (size_t i = 0; i < stoichiometric.size(); i++)
        {
            const auto& reaction = stoichiometric.at(i);
            const auto& indexes = reaction.at(0);
            const auto& nu_fwd = reaction.at(1);
            const auto& nu_bwd = reaction.at(2);
            const auto& nu_vals = reaction.at(3);

            stoich_t sum_nu_reaction = 0;
            for (size_t j = 0; j < num_max_species_per_reaction; j++)
            {
                index_t index = indexes.at(j).get<index_t>();
                stoich_t nu_f = nu_fwd.at(j).get<stoich_t>();
                stoich_t nu_b = nu_bwd.at(j).get<stoich_t>();
                stoich_t nu_val = nu_vals.at(j).get<stoich_t>();

                nu_indexes_temp[i * num_max_species_per_reaction + j]   = index;
                nu_forward_temp[i * num_max_species_per_reaction + j]  = nu_f;
                nu_backward_temp[i * num_max_species_per_reaction + j] = nu_b;
                nu_temp[i * num_max_species_per_reaction + j]          = nu_val;

                sum_nu_reaction += nu_val;
            }
            sum_nu_temp[i] = sum_nu_reaction;
        }
        nu_indexes   = (const index_t*)nu_indexes_temp;
        nu_forward  = (const stoich_t*)nu_forward_temp;
        nu_backward = (const stoich_t*)nu_backward_temp;
        nu          = (const stoich_t*)nu_temp;
        sum_nu      = (const stoich_t*)sum_nu_temp;


        // Falloff reaction types
        if (!check_key_exists(j, "falloff_reactions", "falloff_reaction_types")) return;
        const nlohmann::json& falloff_types = j.at("falloff_reactions").at("falloff_reaction_types");
        if (falloff_types.size() != num_falloff_reactions) return;
        
        reac_type* falloff_reac_types_temp = new reac_type[num_falloff_reactions];
        for (size_t i = 0; i < falloff_types.size(); i++)
        {
            const auto& reaction = falloff_types.at(i);
            std::string type_str = reaction.get<std::string>();
            if (type_str == "lindemann")
                falloff_reac_types_temp[i] = reac_type::lindemann;
            else if (type_str == "troe")
                falloff_reac_types_temp[i] = reac_type::troe;
            else if (type_str == "sri")
                falloff_reac_types_temp[i] = reac_type::sri;
            else
            {
                LOG_ERROR(
                    Error::severity::warning, Error::type::preprocess,
                    "Unknown falloff reaction type \"" + type_str + "\" for reaction index " + std::to_string(i) + ". Using Lindemann value.");
                falloff_reac_types_temp[i] = reac_type::lindemann;
            }
        }
        falloff_reaction_types = (const reac_type*)falloff_reac_types_temp;

    }
    catch(const std::exception& e)
    {
        LOG_ERROR(
            Error::severity::error,
            Error::type::preprocess,
            "Exception occurred while parsing JSON data: " + std::string(e.what())
        );
    }
}


// Removes "//" comments from a multiline string
// Does not support the placement of "//" within the JSON itself
std::string strip_comments(const std::string& s) {
    std::string out;
    out.reserve(s.size());  // pre-allocate to avoid reallocations
    
    const size_t n = s.size();
    for (size_t i = 0; i < n; ++i) {
        char c = s[i];
        
        // Check for "//" comment
        if (c == '/' && i + 1 < n && s[i + 1] == '/') {
            // Skip until end of line
            i += 2;
            while (i < n && s[i] != '\n' && s[i] != '\r') {
                ++i;
            }
            // Keep the newline character (preserves line numbers for error messages)
            if (i < n) {
                out.push_back('\n');
            }
        } else {
            out.push_back(c);
        }
    }

    return out;
}


json get_json(const std::string& json_path)
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

    // Read text from file, strip comments
    std::string raw_content(
        (std::istreambuf_iterator<char>(input_file)),
        std::istreambuf_iterator<char>()
    );
    input_file.close();
    std::string cleaned = strip_comments(raw_content);

    // Read JSON data
    nlohmann::json j;
    try
    {
        j = json::parse(cleaned);
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

    // Validate JSON content
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


Parameters::~Parameters()
{
    if (W != nullptr) delete[] W;
    if (thermal_conductivities != nullptr) delete[] thermal_conductivities;
    if (sqrt_van_der_waals_a != nullptr) delete[] sqrt_van_der_waals_a;
    if (van_der_waals_b != nullptr) delete[] van_der_waals_b;
    if (temp_ranges != nullptr) delete[] temp_ranges;
    if (a_low != nullptr) delete[] a_low;
    if (a_high != nullptr) delete[] a_high;
    if (interval_values != nullptr) delete[] interval_values;
    if (interval_derivatives != nullptr) delete[] interval_derivatives;
    if (arrhenius_parameters != nullptr) delete[] arrhenius_parameters;
    if (nu_indexes != nullptr) delete[] nu_indexes;
    if (nu_forward != nullptr) delete[] nu_forward;
    if (nu_backward != nullptr) delete[] nu_backward;
    if (reaction_order != nullptr) delete[] reaction_order;
    if (nu != nullptr) delete[] nu;
    if (sum_nu != nullptr) delete[] sum_nu;
    if (third_body_reaction_indexes != nullptr) delete[] third_body_reaction_indexes;
    if (is_falloff_reaction != nullptr) delete[] is_falloff_reaction;
    if (alpha != nullptr) delete[] alpha;
    if (irreversible_reaction_indexes != nullptr) delete[] irreversible_reaction_indexes;
    if (falloff_reaction_indexes != nullptr) delete[] falloff_reaction_indexes;
    if (falloff_reaction_types != nullptr) delete[] falloff_reaction_types;
    if (is_third_body_indexes != nullptr) delete[] is_third_body_indexes;
    if (falloff_parameters != nullptr) delete[] falloff_parameters;
    if (troe_parameters != nullptr) delete[] troe_parameters;
    if (sri_parameters != nullptr) delete[] sri_parameters;
    if (plog_reaction_indexes != nullptr) delete[] plog_reaction_indexes;
    if (plog_seperators != nullptr) delete[] plog_seperators;
    if (plog_parameters != nullptr) delete[] plog_parameters;
}


const Parameters *Parameters::get_parameters(const std::string& mech_name)
{
    // Check if already loaded (by Parameters::mechanism_name)
    Timer timer; timer.start();
    auto it = std::find_if(_mechanisms.begin(), _mechanisms.end(),
        [&](const std::unique_ptr<Parameters>& p) {
            return p && p->mechanism_name == mech_name;
        });
    if (it != _mechanisms.end())
        return it->get();
    
    // Check if JSON directory exists
    std::filesystem::path mech_dir_path = std::filesystem::path("mechanism") / "json_files";
    if (!std::filesystem::exists(mech_dir_path))
    {
        LOG_ERROR(
            Error::severity::error, Error::type::preprocess,
            "Directory \"" + mech_dir_path.string() + "\" not found. Make sure, your current working directory is correct."
        );
        return nullptr;
    }

    // Check if corresponding JSON file exists
    std::filesystem::path json_path = mech_dir_path / (mech_name + ".json");
    if (!std::filesystem::exists(json_path))
    {
        // Collect available mechanisms (JSON files without extension)
        std::vector<std::string> available;
        for (const auto& entry : std::filesystem::directory_iterator(mech_dir_path))
        {
            if (!entry.is_regular_file()) continue;
            const auto& p = entry.path();
            if (p.has_extension() && p.extension() == ".json")
                available.emplace_back(p.stem().string());
        }

        LOG_ERROR(
            Error::severity::error, Error::type::preprocess,
            "Mechanism JSON file \"" + json_path.string() + "\" not found. Check if the mechanism name is correct. Available: " + ::to_string((std::string*)available.data(), available.size())
        );
        return nullptr;
    }
    
    // Load new Parameters instance
    auto param = std::make_unique<Parameters>(json_path.string());
    const double runtime = timer.lap();
    LOG_ERROR(
        Error::severity::info, Error::type::preprocess,
        "Loaded mechanism \"" + param->mechanism_name + "\" from \"" + json_path.string() + "\" in " + Timer::format_time(runtime) + "."
    );
    const std::string expected_mechanism_name = json_path.stem().string();
    if (param->mechanism_name != expected_mechanism_name)
    {
        LOG_ERROR(
            Error::severity::error, Error::type::preprocess,
            "Mechanism name in JSON file \"" + param->mechanism_name + "\" does not match the expected name from the file name \"" + expected_mechanism_name + "\"."
        );
        return nullptr;
    } 

    _mechanisms.push_back(std::move(param));
    return _mechanisms.back().get();
}


index_t Parameters::get_species(std::string name) const
{
    std::transform(name.begin(), name.end(), name.begin(), ::toupper);
    auto it = _species.find(name);
    if (it == _species.end())
    {
        std::stringstream ss;
        ss << "Species \"" << name << "\" not found in " << this->mechanism_name << ". Valid species are: ";
        ss << ::to_string((std::string*)Parameters::species_names.data(), Parameters::species_names.size());
        LOG_ERROR(Error::severity::error, Error::type::preprocess, ss.str());
        return invalid_index;
    }
    return it->second;
}


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
    return Parameters::excitation::sinusoid; // Default to a known excitation type
}
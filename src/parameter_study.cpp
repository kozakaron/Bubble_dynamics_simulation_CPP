#include <iomanip>
#include <sstream>
#include <numbers>
#include <cmath>
#include <regex>
#include <limits>
#include <filesystem>

#include "nlohmann/json.hpp"
#include "parameter_study.h"
#include "ode_fun.h"
#include "ode_solver_sundials.h"

using ordered_json = nlohmann::ordered_json;


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

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
#include "ode_solver.h"




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

    return save_folder_path.string();
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
        nlohmann::ordered_json j = parameter_combinator.to_json();
        j["info"] = nlohmann::ordered_json::object({
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
    OdeSolver solver(4 + par->num_species);

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
        SimulationData data = solver.solve(
            this->t_max,                // t_max [s]
            &ode,                       // ode_ptr
            this->timeout,              // timeout [s]
            false                       // save solution
        );

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

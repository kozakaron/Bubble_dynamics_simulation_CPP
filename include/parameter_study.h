#ifndef PARAMETER_STUDY_H
#define PARAMETER_STUDY_H
#include "nlohmann/json_fwd.hpp"
#include "common.h"
#include "parameters.h"
#include "control_parameters.h"
#include "ode_solver.h"
#include "parameter_combinator.h"

#include <string>
#include <ostream>
#include <fstream>
#include <atomic>
#include <mutex>
#include <thread>

class ParameterStudy
{
private:
    ParameterCombinator &parameter_combinator;
    std::string save_folder;
    std::ofstream output_log_file;
    std::mutex output_mutex;
    double best_energy_demand;
    const double t_max;
    const double timeout;
    std::atomic<size_t> successful_simulations;
    std::atomic<size_t> total_simulations;
    
    void parameter_study_task(const bool print_output, const size_t thread_id);
public:
    ParameterStudy(
        ParameterCombinator &parameter_combinator,
        std::string save_folder,
        const double t_max = 1.0,
        const double timeout = 60.0
    );
    ~ParameterStudy();
    void run(const size_t num_threads, const bool print_output=true);
};

#endif // PARAMETER_STUDY_H
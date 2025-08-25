#ifndef ODE_SOLUTION_H
#define ODE_SOLUTION_H

#include "nlohmann/json_fwd.hpp"
#include "common.h"
#include "control_parameters.h"
#include "parameter_combinator.h"

#include <vector>
#include <string>
#include <ostream>
#include <fstream>


class OdeSolution
{
public:
// Static members
    static constexpr char csv_header[] = "success,num_dim,num_steps,num_repeats,num_fun_evals,num_fun_evals_jac,num_jac_evals,num_lin_iters,num_nonlin_iters,runtime,t_last,x_0,x_last";
// Solution and statistics
    std::vector<double> t;
    std::vector<std::vector<double>> x;
    size_t num_dim;
    size_t num_steps;
    size_t num_repeats;
    size_t num_fun_evals;
    size_t num_fun_evals_jac;
    size_t num_jac_evals;
    size_t num_lin_iters;
    size_t num_nonlin_iters;
    std::vector<double> total_error;
// Metadata
    double runtime;
    size_t error_ID;
// Methods
    OdeSolution();
    ~OdeSolution();
    is_success success() const;
    void push_t_x(const double t_i, const double *x_i);
    void clear();
    std::string to_csv() const;
    std::string to_string(const bool colored=true, const bool with_code=true) const;
    nlohmann::ordered_json to_json() const;
    friend std::ostream &operator<<(std::ostream &os, const OdeSolution &ode);
};


class SimulationData
{
public:
// Static members
    static const std::string csv_header;
    static const Error no_error;
    static const double infinite_energy_demand;
// Mid-processing results
    double R_max;              // [m]
    double T_max;              // [K]
    double t_peak;             // [s]
    double R_min;              // [m]
    double T_min;              // [K]
    double t_collapse;         // [s]
// Post-processing results
    double dissipated_energy;  // [J]
    double n_target_specie;    // [mol]
    double energy_demand;      // [MJ/kg]
// Members
    const ControlParameters &cpar;
    OdeSolution sol;
// Methods
    SimulationData(const ControlParameters &cpar);
    void midprocess(const double t, const double* x);
    void postprocess();
    std::string to_csv() const;
    std::string to_string() const;
    std::string to_small_string(const ParameterCombinator &ps, const double best_energy_demand, const bool colored=true) const;
    nlohmann::ordered_json to_json() const;
    void save_json_with_binary(const std::string &json_path) const;
    friend std::ostream &operator<<(std::ostream &os, const SimulationData &data);
};


#endif // ODE_SOLUTION_H
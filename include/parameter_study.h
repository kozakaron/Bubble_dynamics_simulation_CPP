#ifndef PARAMETER_STUDY_H
#define PARAMETER_STUDY_H
#include "common.h"
#include "parameters.h"
#include "control_parameters.h"
#include "ode_solver.h"

#include <string>
#include <ostream>
#include <fstream>
#include <memory>
#include <vector>
#include <variant>
#include <atomic>
#include <mutex>
#include <functional>

// Represents possible values for a parameter in a parameter study. Options:
// - Const: constant value
// - LinearRange: linear subdivision of an interval
// - PowRange: uneven subdivision of an interval
class Range
{
protected:
    double start;
    double end;
    size_t num_steps;
public:
    Range();
    Range(double start, double end, size_t num_steps);
    virtual ~Range() = default;
    virtual std::string to_string() const = 0;
    virtual std::string to_array() const;
    friend std::ostream &operator<<(std::ostream &os, const Range &range);
    virtual double operator[](size_t i) const = 0;
    size_t get_num_steps() const;
};

// Represents a constant value. f(x) = c
// E.g.: .alfa =  Const(0.35);   // alfa = {0.35}
class Const: public Range
{
public:
    Const(double value);
    double operator[](size_t i) const override;
    std::string to_string() const override;
};

// Represents an even subdivision of an interval. f(x) = a * x + b
// E.g.: .R_E = LinearRange(1e-6, 10e-6, 5);   // R_E = {1e-6, 3.25e-6, 5.5e-6, 7.75e-6, 1e-5}
class LinearRange : public Range
{
public:
    LinearRange(double start, double end, size_t num_steps);
    double operator[](size_t i) const override;
    std::string to_string() const override;
};

// Represents an uneven subdivision of an interval. f(x) = a * x^b + c
// If 0 < base < 1, steps are getting gradually smaller.
// If 1 < base, steps are getting larger. Increase base to get a larger increase in steps.
// E.g.: .R_E = PowRange(1e-6, 10e-6, 5, 2);   // R_E = {1.0000e-06, 1.5625e-06, 3.2500e-06, 6.0625e-06, 1.0000e-05}
class PowRange : public Range
{
private:
    double base;
    double a;
    double c;
public:
    PowRange(double start, double end, size_t num_steps, double base=2.0);
    double operator[](size_t i) const override;
    std::string to_string() const override;
};


class ParameterCombinator
{
private:
    std::atomic<size_t> combination_ID;
    size_t total_combination_count;

    Parameters::mechanism mechanism;
    std::unique_ptr<Range> R_E;
    std::vector<std::string> species;
    std::vector<double> fractions;
    std::unique_ptr<Range> P_amb;
    std::unique_ptr<Range> T_inf;
    std::unique_ptr<Range> alfa_M;
    std::unique_ptr<Range> P_v;
    std::unique_ptr<Range> mu_L;
    std::unique_ptr<Range> rho_L;
    std::unique_ptr<Range> c_L;
    std::unique_ptr<Range> surfactant;
    bool enable_heat_transfer;
    bool enable_evaporation;
    bool enable_reactions;
    bool enable_dissipated_energy;
    std::string target_specie;
    std::vector<std::unique_ptr<Range>> excitation_params;
    Parameters::excitation excitation_type;
public:
    typedef std::variant<Const, LinearRange, PowRange> AnyRange;
    struct Builder {
        Parameters::mechanism mechanism         = Parameters::mechanism::chemkin_ar_he;
        AnyRange R_E                            = Const(10.0e-6);
        std::vector<std::string> species        = {"O2"};
        std::vector<double> fractions           = {1.0};
        AnyRange P_amb                          = Const(101325.0);
        AnyRange T_inf                          = Const(293.15);
        AnyRange alfa_M                         = Const(0.35);
        AnyRange P_v                            = Const(2338.1);
        AnyRange mu_L                           = Const(0.001);
        AnyRange rho_L                          = Const(998.2);
        AnyRange c_L                            = Const(1483.0);
        AnyRange surfactant                     = Const(1.0);
        bool enable_heat_transfer               = true;
        bool enable_evaporation                 = true;
        bool enable_reactions                   = true;
        bool enable_dissipated_energy           = true;
        std::string target_specie               = "H2";
        std::vector<AnyRange> excitation_params = {Const(-2.0e5), Const(30000.0), Const(1.0)};
        Parameters::excitation excitation_type  = Parameters::excitation::sin_impulse;
    };

    ParameterCombinator(const Builder &builder);
    std::string to_string(const bool with_code=false) const;
    friend std::ostream &operator<<(std::ostream &os, const ParameterCombinator &pc);
    friend class SimulationData;
    friend class ParameterStudy;
    size_t get_total_combination_count() const;
    size_t get_next_combination_ID() const;
    std::pair<is_success, ControlParameters> get_next_combination();
    const Parameters* get_mechanism_parameters() const;
};


class SimulationData
{
public:
    static const std::string csv_header;
    static const Error no_error;
    static const double infinite_energy_demand;

    const ControlParameters &cpar;
    const OdeSolution &sol;
    const double dissipated_energy;  // [J]
    const double n_target_specie;    // [mol]
    const double energy_demand;      // [MJ/kg]

    SimulationData(
        const ControlParameters &cpar,
        const OdeSolution &sol
    );
    std::string to_csv() const;
    std::string to_string() const;
    std::string to_small_string(const ParameterCombinator &ps, const double best_energy_demand, const bool colored=true) const;
    friend std::ostream &operator<<(std::ostream &os, const SimulationData &data);
};


class ParameterStudy
{
private:
    ParameterCombinator &parameter_combinator;
    std::string save_folder;
    std::ofstream output_log_file;
    std::mutex output_mutex;
    std::function<OdeSolver*(size_t)> solver_factory;
    double best_energy_demand;
    const double t_max;
    const double timeout;
    
    void parameter_study_task(const bool print_output, const size_t thread_id);
public:
    ParameterStudy(
        ParameterCombinator &parameter_combinator,
        std::string save_folder,
        std::function<OdeSolver*(size_t)> solver_factory,
        const double t_max = 1.0,
        const double timeout = 60.0
    );
    ~ParameterStudy();
    void run(const size_t num_threads, const bool print_output=true);
};

#endif // PARAMETER_STUDY_H
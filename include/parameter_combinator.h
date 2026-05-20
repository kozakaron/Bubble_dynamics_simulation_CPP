#ifndef PARAMETER_COMBINATOR_H
#define PARAMETER_COMBINATOR_H
#include "nlohmann/json_fwd.hpp"
#include "common.h"
#include "parameters.h"
#include "control_parameters.h"

#include <string>
#include <ostream>
#include <memory>
#include <vector>
#include <variant>
#include <atomic>

// Represents possible values for a parameter in a parameter study. Options:
// - Const: constant value
// - LinearRange: linear subdivision of an interval
// - LogRange: uneven subdivision of an interval, same as numpy.logspace in Python
// - GeomRange: uneven subdivision of an interval, more versatile than LogRange
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
    virtual nlohmann::ordered_json to_json() const = 0;
    virtual std::string to_array() const;
    friend std::ostream &operator<<(std::ostream &os, const Range &range);
    virtual double operator[](size_t i) const = 0;
    size_t get_num_steps() const;
};

// Represents a constant value. f(x) = c
// E.g.: .alpha =  Const(0.35);   // alpha = {0.35}
class Const: public Range
{
public:
    Const(double value);
    double operator[](size_t i) const override;
    std::string to_string() const override;
    nlohmann::ordered_json to_json() const override;
};

// Represents an even subdivision of an interval. f(x) = a * x + b
// E.g.: .R_E = LinearRange(1e-6, 10e-6, 5);   // R_E = {1e-6, 3.25e-6, 5.5e-6, 7.75e-6, 1e-5}
class LinearRange : public Range
{
public:
    LinearRange(double start, double end, size_t num_steps);
    double operator[](size_t i) const override;
    std::string to_string() const override;
    nlohmann::ordered_json to_json() const override;
};

// Represents an uneven subdivision of an interval. f(x) = exp(a + b * x)
// Same as numpy.logspace in Python.
// E.g.: .R_E = LogRange(1e-6, 10e-6, 5);   // R_E = {1e-06, 1.77828e-06, 3.16228e-06, 5.62341e-06, 1e-05}
class LogRange : public Range
{
private:
    double a;
    double b;
    int sign;
public:
    LogRange(double start, double end, size_t num_steps);
    double operator[](size_t i) const override;
    std::string to_string() const override;
    nlohmann::ordered_json to_json() const override;
};

// Represents an uneven subdivision of an interval. f(x) = start + a_1 + a_1*q + a_1*q^2 + ... + a_1*q^(x-1) = start + a_1 * (1 - q^x) / (1 - q)
// If 0 < q < 1, steps are getting gradually smaller.
// If 1 < q, steps are getting larger. Increase q to get a larger increase in steps.
// E.g.: .R_E = GeomRange(1e-6, 10e-6, 5, 2.0);   // R_E = {1e-06, 1.6e-06, 2.8e-06, 5.2e-06, 1e-05}
class GeomRange : public Range
{
private:
    double a_1;
    double q;
public:
    GeomRange(double start, double end, size_t num_steps, double q);
    double operator[](size_t i) const override;
    std::string to_string() const override;
    nlohmann::ordered_json to_json() const override;
};

class ParameterCombinator
{
private:
    std::atomic<size_t> combination_ID;
    size_t total_combination_count;

    std::string mechanism;
    std::unique_ptr<Range> R_E;
    std::unique_ptr<Range> ratio;
    std::vector<std::string> species;
    std::vector<double> fractions;
    std::unique_ptr<Range> P_amb;
    std::unique_ptr<Range> T_inf;
    std::unique_ptr<Range> alpha_M;
    std::unique_ptr<Range> P_v;
    std::unique_ptr<Range> mu_L;
	std::unique_ptr<Range> nu_L;
	std::unique_ptr<Range> Gamma_L;
	std::unique_ptr<Range> B_L;
	std::unique_ptr<Range> b_L;
    std::unique_ptr<Range> rho_0;
    std::unique_ptr<Range> c_0;
	std::unique_ptr<Range> sigma_var;
    std::unique_ptr<Range> surfactant;
    bool enable_heat_transfer;
    bool enable_evaporation;
    bool enable_reactions;
    bool enable_dissipated_energy;
    bool enable_van_der_waals;
    bool enable_rate_thresholding;
    std::string target_specie;
    Parameters::excitation excitation_type;
    std::vector<std::unique_ptr<Range>> excitation_params;
    std::unique_ptr<Range> excitation_cycles;
    std::unique_ptr<Range> ramp_up_cycles;
	std::unique_ptr<Range> R_and_R_dot_from_file;
	std::unique_ptr<Range> rows;
	std::unique_ptr<Range> cols;
	std::vector<std::string> file_name;
	std::unique_ptr<Range> const_V;
	std::unique_ptr<Range> const_T;
	std::vector<std::string> EoS_liquid;
	
public:
    typedef std::variant<Const, LinearRange, LogRange, GeomRange> AnyRange;
    struct Builder {
        std::string mechanism                   = "chemkin_elte2016_hydrogen";
        AnyRange R_E                            = Const(10.0e-6);
        AnyRange ratio                          = Const(1.0);
        std::vector<std::string> species        = {"O2"};
        std::vector<double> fractions           = {1.0};
        AnyRange P_amb                          = Const(101325.0);
        AnyRange T_inf                          = Const(293.15);
        AnyRange alpha_M                        = Const(0.35);
        AnyRange P_v                            = Const(2338.1);
        AnyRange mu_L                           = Const(0.001);
		AnyRange nu_L                           = Const(1.0);
		AnyRange Gamma_L                        = Const(8.46);
		AnyRange B_L                            = Const(616640000.0);
		AnyRange b_L                            = Const(6.72e-4);
        AnyRange rho_0                          = Const(998.2);
        AnyRange c_0                            = Const(1483.0);
		AnyRange sigma_var                      = Const(0.0728);
        AnyRange surfactant                     = Const(1.0);
        bool enable_heat_transfer               = true;
        bool enable_evaporation                 = true;
        bool enable_reactions                   = true;
        bool enable_dissipated_energy           = true;
        bool enable_van_der_waals               = true;
        bool enable_rate_thresholding           = true;
        std::string target_specie               = "H2";
        Parameters::excitation excitation_type  = Parameters::excitation::sinusoid;
        std::vector<AnyRange> excitation_params = {Const(-2.0e5), Const(30000.0)};	
        AnyRange excitation_cycles              = Const(1.0);
        AnyRange ramp_up_cycles                 = Const(0.0);
		AnyRange R_and_R_dot_from_file 			= Const(0);
		AnyRange rows 							= Const(1);
		AnyRange cols 							= Const(1);
		std::string file_name 					= "temp.csv";
		AnyRange const_V						= Const(0);
		AnyRange const_T						= Const(0);
		std::string EoS_liquid 					= "Tait";
    };

    ParameterCombinator(const Builder &builder);
    ParameterCombinator(const nlohmann::ordered_json& j);
    ParameterCombinator(const std::string& json_path);
    std::string to_string(const bool with_code=false) const;
    nlohmann::ordered_json to_json() const;
    friend std::ostream &operator<<(std::ostream &os, const ParameterCombinator &pc);
    friend class SimulationData;
    friend class ParameterStudy;
    size_t get_total_combination_count() const;
    size_t get_next_combination_ID() const;
    std::pair<is_success, ControlParameters> get_next_combination();
    const Parameters* get_mechanism_parameters() const;

private:
    void init(const ParameterCombinator::Builder &builder);
    void init(const nlohmann::ordered_json &j);
};


#endif // PARAMETER_COMBINATOR_H
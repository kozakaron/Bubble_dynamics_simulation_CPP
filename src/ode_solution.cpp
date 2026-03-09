#include <sstream>
#include <iomanip>
#include <numbers>
#include <cmath>
#include <regex>
#include <limits>
#include <filesystem>

#include "nlohmann/json.hpp"
#include "ode_solution.h"
#include "ode_fun.h"


OdeSolution::OdeSolution():
    t(),
    x(),
    num_dim(0),
    num_steps(0),
    num_repeats(0),
    num_fun_evals(0),
    num_fun_evals_jac(0),
    num_jac_evals(0),
    num_lin_iters(0),
    num_nonlin_iters(0),
    total_error(),
    runtime(0.0),
    error_ID(ErrorHandler::no_error)
{

}

OdeSolution::~OdeSolution() { }

is_success OdeSolution::success() const
{
    return error_ID == ErrorHandler::no_error;
}


void OdeSolution::push_t_x(const double t_dimless, const double *x_dimless, OdeFun* ode_ptr)
{
    t.push_back(t_dimless);
    x.push_back(std::vector<double>(x_dimless, x_dimless + num_dim));
    ode_ptr->cpar.dimensionalize(t.back(), x.back().data());
    
    // Compute sum of molar concentrations
    const double T = x.back()[2];
    double M = 0.0;
    for (size_t k = 3; k < num_dim; ++k) {
        M += x.back()[k];
    }
    const double* conc = x.back().data() + 3;
    
    // Compute internal pressure (using dimensional values)
    double p_int = ode_ptr->internal_pressure(T, M, conc);
    p_internal.push_back(p_int);
    
    // Compute excitation pressure (external pressure)
    auto [P_inf, P_inf_dot] = ode_ptr->excitation_pressures(t.back());
    (void)P_inf_dot;
    p_excitation.push_back(P_inf);
}

void OdeSolution::clear()
{
    t.clear();
    x.clear();
    p_excitation.clear();
    p_internal.clear();
    num_dim = 0;
    num_steps = 0;
    num_repeats = 0;
    num_fun_evals = 0;
    num_fun_evals_jac = 0;
    num_jac_evals = 0;
    num_lin_iters = 0;
    num_nonlin_iters = 0;
    total_error.clear();
    runtime = 0.0;
    error_ID = ErrorHandler::no_error;
}


std::string percent(const size_t num, const size_t num_steps)
{
    std::stringstream ss;
    ss << "    // ";
    if (num_steps == 0)
    {
        ss << std::setw(8) << "nan %";
    }
    else if (num == num_steps)
    {
        ss << std::setw(8) << 100.0 << " %";
    }
    else if (num < num_steps)
    {
        ss << std::setw(8) << std::setprecision(2) << 100.0 * num / num_steps << " %";
    } else
    {
        ss << std::setw(8) << std::setprecision(2) << (double) num / num_steps << " per step";
    }
    
    return ss.str();
}


std::string OdeSolution::to_csv() const
{
    std::stringstream ss;
    auto format_double = [](std::ostream& os) -> std::ostream& {
        return os << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
    };

    ss << (this->success() ? "true" : "false") << "," << this->num_dim << "," << this->num_steps << ",";
    ss << this->num_repeats << "," << this->num_fun_evals << "," << this->num_fun_evals_jac << "," << this->num_jac_evals << ",";
    ss << this->num_lin_iters << "," << this->num_nonlin_iters << ",";
    ss << format_double << this->runtime << "," << format_double;
    if (this->t.empty()) ss << format_double << 0.0 << ",";
    else ss << this->t.back() << ",";
    if (!this->x.empty())
    {
        for (size_t i = 0; i < this->x[0].size(); ++i)
            ss << format_double << this->x[0][i] << ";";
        ss << ",";
        for (size_t i = 0; i < this->x.back().size(); ++i)
            ss << format_double << this->x.back()[i] << ";";
    }
    else
    {
        ss << ",";
    }
    
    return ss.str();
}


std::string OdeSolution::to_string(const bool colored, const bool with_code) const
{
    std::stringstream ss;
    const size_t strw = 24;
    const size_t intw = 10;
    
    ss << std::left;
    if (with_code) ss << "OdeSolution{\n";
    ss << std::setw(strw) << "    .success" << " = ";
    if (success())
    {
        if (colored) ss << colors::green << colors::bold;
        ss << "true";
        if (colored) ss << colors::reset;
        ss << ",\n";
    }
    else
    {
        if (colored) ss << colors::red << colors::bold;
        ss << "false";
        if (colored) ss << colors::reset;
        ss << ",\n";
        ss << std::setw(strw) << "    .error" << " = \"" << ErrorHandler::get_error(error_ID).to_string(colored) << "\",\n";
    }

    ss << std::setw(strw) << "    .runtime"            << " = \"" << Timer::format_time(runtime) << "\",\n";
    ss << std::setw(strw) << "    .num_steps"          << " = " << std::setw(intw) << std::to_string(num_steps)          + "," << percent(num_steps, num_steps) << "\n";
    ss << std::setw(strw) << "    .num_saved_steps"    << " = " << std::setw(intw) << std::to_string(x.size())           + "," << percent(x.size(), num_steps) << "\n";
    ss << std::setw(strw) << "    .num_repeats"        << " = " << std::setw(intw) << std::to_string(num_repeats)        + "," << percent(num_repeats, num_steps) << "\n";
    ss << std::setw(strw) << "    .num_fun_evals"      << " = " << std::setw(intw) << std::to_string(num_fun_evals)      + "," << percent(num_fun_evals, num_steps) << "\n";
    ss << std::setw(strw) << "    .num_fun_evals_jac"  << " = " << std::setw(intw) << std::to_string(num_fun_evals_jac)  + "," << percent(num_fun_evals_jac, num_steps) << "\n";
    ss << std::setw(strw) << "    .num_jac_evals"      << " = " << std::setw(intw) << std::to_string(num_jac_evals)      + "," << percent(num_jac_evals, num_steps) << "\n";
    ss << std::setw(strw) << "    .num_lin_iters"      << " = " << std::setw(intw) << std::to_string(num_lin_iters)      + "," << percent(num_lin_iters, num_steps) << "\n";
    ss << std::setw(strw) << "    .num_nonlin_iters"   << " = " << std::setw(intw) << std::to_string(num_nonlin_iters)   + "," << percent(num_nonlin_iters, num_steps) << "\n";
    if (this->t.size() >=2 && this->x.size() >= 2 && this->num_dim > 0)
    {
        ss << std::setw(strw) << "    .t" << " = " << "{";
        ss << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10) << this->t[0]                  << ", /* ..., */ ";
        ss << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10) << this->t[this->t.size() - 1] << "},\n";
        ss << std::setw(strw) << "    .x" << " = " << "{\n";
        ss << std::setw(strw+8) << " " << ::to_string((double*)this->x[0].data(),                  this->x[0].size())                  << ",\n";
        ss << std::setw(strw+8) << " " << " // ..., \n";
        ss << std::setw(strw+8) << " " << ::to_string((double*)this->x[this->x.size() - 1].data(), this->x[this->x.size() - 1].size()) << "\n";
        ss << std::setw(strw+4) << " " << "},\n";
    }
    ss << std::setw(strw) << "    .total_error"        << " = " << "{" << ::to_string((double*)this->total_error.data(), this->total_error.size()) << "}\n";

    if (with_code) ss << "}";
    ss << std::right;
    return ss.str();
}


nlohmann::ordered_json OdeSolution::to_json() const
{
    nlohmann::ordered_json j;
    j["success"] = this->success();
    j["error"] = ErrorHandler::get_error(this->error_ID).to_string();
    j["runtime"] = this->runtime;
    j["num_dim"] = this->num_dim;
    j["num_steps"] = this->num_steps;
    j["num_saved_steps"] = this->x.size();
    j["num_repeats"] = this->num_repeats;
    j["num_fun_evals"] = this->num_fun_evals;
    j["num_fun_evals_jac"] = this->num_fun_evals_jac;
    j["num_jac_evals"] = this->num_jac_evals;
    j["num_lin_iters"] = this->num_lin_iters;
    j["num_nonlin_iters"] = this->num_nonlin_iters;
    j["total_error"] = this->total_error;
    // Saved as binary data in SimulationData:
    //j["t"] = std::vector<double>({this->t.front(), this->t.back()});
    //j["x"] = std::vector<std::vector<double>>({this->x.front(), this->x.back()});
    //j["p_excitation"] = std::vector<double>({this->p_excitation.front(), this->p_excitation.back()});
    //j["p_internal"] = std::vector<double>({this->p_internal.front(), this->p_internal.back()});
    
    return j;
}



std::ostream &operator<<(std::ostream &os, const OdeSolution &ode)
{
    os << ode.to_string();
    return os;
}


const std::string SimulationData::csv_header = std::string("dissipated_energy,expansion_work,n_target_specie,energy_demand,R_max,T_max,t_peak,R_min,T_min,v_max,p_internal_max,p_internal_min,Ma_max,T_L_max,c_L_max,rho_L_max,")
                                             + std::string(ControlParameters::csv_header) + std::string(",")
                                             + std::string(OdeSolution::csv_header) + std::string(",") + std::string(Error::csv_header);

const Error SimulationData::no_error = Error(Error::severity::info, Error::type::general, "No error", "", __FILE__, __LINE__, 0);

const double SimulationData::infinite_energy_demand = std::numeric_limits<double>::max();

SimulationData::SimulationData(const ControlParameters &cpar):
    R_max(0.0),
    T_max(0.0),
    t_peak(0.0),
    R_min(std::numeric_limits<double>::infinity()),
    T_min(std::numeric_limits<double>::infinity()),
    v_max(0.0),
    p_internal_max(0.0),
    p_internal_min(cpar.P_amb),
    Ma_max(0.0),
    T_L_max(cpar.T_inf),
    c_L_max(cpar.c_L),
    rho_L_max(cpar.rho_L),
    dissipated_energy(0.0),
    expansion_work(0.0),
    n_target_specie(0.0),
    energy_demand(std::numeric_limits<double>::infinity()),
    cpar(cpar),
    sol(),
    x_dimensional()
{}


void SimulationData::midprocess(const double t_dimless, const double* x_dimless, OdeFun* ode_ptr)
{
    const Parameters* par = cpar.par;
    const size_t num_dim = (par != nullptr) ? (4 + par->num_species) : 4;

    // Dimensionalize a local copy of the state
    double t_dimensional = t_dimless;
    x_dimensional.assign(x_dimless, x_dimless + num_dim);
    cpar.dimensionalize(t_dimensional, x_dimensional.data());

    const double R     = x_dimensional[0];       // [m]
    const double R_dot = x_dimensional[1];       // [m/s]
    const double T     = x_dimensional[2];       // [K]
    const double* conc = x_dimensional.data() + 3;
    double M = 0.0;
    for (size_t k = 3; k < num_dim; ++k) {
        M += x_dimensional[k];
    }
    const double p_int = ode_ptr->internal_pressure(T, M, conc);
    auto [P_inf, P_inf_dot] = ode_ptr->excitation_pressures(t_dimensional);
    const double p_L = p_int - (2.0 * cpar.surfactant * par->sigma + 4.0 * cpar.mu_L * R_dot) / R;
    auto [c_L, rho_L, rho_inf, H_unused] = ode_ptr->liquid_properties(p_L, P_inf);
    const double Ma = std::abs(R_dot) / c_L;
    constexpr double Gamma_L = Parameters::nasg::Gamma_L;
    constexpr double B_L = Parameters::nasg::B_L;
    const double T_L = cpar.T_inf * std::pow((p_L + B_L) / (cpar.P_amb + B_L), (Gamma_L - 1.0) / Gamma_L); // Only for NASG
    (void)P_inf_dot;
    (void)rho_inf; (void)H_unused;

    // Update fields
    if (R < R_min) R_min = R;
    if (R > R_max) R_max = R;
    if (T < T_min) T_min = T;
    if (T > T_max)
    {
        T_max = T;
        t_peak = t_dimensional;
    }
    if (std::abs(R_dot) > v_max) v_max = std::abs(R_dot);
    if (p_int > p_internal_max) p_internal_max = p_int;
    if (p_int < p_internal_min) p_internal_min = p_int;
    if (c_L   > c_L_max)   c_L_max   = c_L;
    if (rho_L > rho_L_max) rho_L_max = rho_L;
    if (Ma > Ma_max) Ma_max = Ma;
    if (T_L > T_L_max) T_L_max = T_L;
}


void SimulationData::postprocess()
{
// dissipated energy [J]
    if (sol.x.empty()) dissipated_energy = 0.0;
    else if (sol.x.back().empty()) dissipated_energy = 0.0;
    else dissipated_energy = sol.x.back().back();

    if (dissipated_energy < -1.0e-8)
        LOG_ERROR(Error::severity::warning, Error::type::postprocess, "Dissipated energy is negative: " + std::to_string(dissipated_energy), cpar.ID);

// expansion work [J]
    const Parameters* par = cpar.par;
    if (
        par == nullptr ||
        cpar.ratio == 1.0
    )
    {
        expansion_work = 0.0;
    }
    else{
        const double R_0 = cpar.ratio * cpar.R_E;   // [m]
        const double V_E = 4.0 / 3.0 * std::numbers::pi * cpar.R_E * cpar.R_E * cpar.R_E;    // [m^3]
        const double V_0 = 4.0 / 3.0 * std::numbers::pi * R_0 * R_0 * R_0;   // [m^3]

        const double p_E = cpar.P_amb + 2.0 * cpar.surfactant * par->sigma / cpar.R_E;   // [Pa]
        const double p_gas = cpar.enable_evaporation ? p_E - cpar.P_v : p_E;   // [Pa]
        const double n_gas = p_gas * V_E / (par->R_g * cpar.T_inf);    // [mol]

        const double W_gas0 = cpar.enable_evaporation ? -(cpar.P_v * V_0 + n_gas * par->R_g * cpar.T_inf * std::log(V_0)) : -(n_gas * par->R_g * cpar.T_inf * std::log(V_0)); // [J]
        const double W_gasE = cpar.enable_evaporation ? -(cpar.P_v * V_E + n_gas * par->R_g * cpar.T_inf * std::log(V_E)) : -(n_gas * par->R_g * cpar.T_inf * std::log(V_E)); // [J]
        const double W_gas = W_gas0 - W_gasE; // [J]

        const double W_surface_tension = par->sigma * cpar.surfactant *  4.0 * std::numbers::pi * (R_0*R_0 - cpar.R_E*cpar.R_E); // [J]
        const double W_flow = cpar.P_amb * 4.0 / 3.0 * std::numbers::pi * (R_0*R_0*R_0 - cpar.R_E*cpar.R_E*cpar.R_E); // [J]
        expansion_work = W_gas + W_surface_tension + W_flow;    // [J]

        if (expansion_work < -1.0e-8)
            LOG_ERROR(Error::severity::warning, Error::type::postprocess, "Expansion work is negative: " + std::to_string(expansion_work), cpar.ID);
    }


// n_target_specie [mol]
    if (
        sol.x.empty() ||
        sol.x.back().empty() ||
        par == nullptr ||
        cpar.target_specie == par->invalid_index ||
        sol.x.front().size() != size_t(4 + par->num_species) ||
        sol.x.back().size() != size_t(4 + par->num_species)
    )
    {
        n_target_specie = 0.0;
    }
    else
    {
        const double R_last = sol.x.back()[0];  // [m]
        const double V_last = 4.0 / 3.0 * std::numbers::pi * std::pow(R_last, 3); // [m^3]
        const double c_target = sol.x.back()[3+cpar.target_specie];  // [mol/m^3]

        n_target_specie = c_target * V_last;  // [mol]
    }

    if (n_target_specie < -1.0e-8)
        LOG_ERROR(Error::severity::warning, Error::type::postprocess, "Target specie concentration is negative: " + std::to_string(n_target_specie), cpar.ID);

// energy demand [MJ/kg]
    const double m_target = n_target_specie * par->W[cpar.target_specie];  // [kg]
    energy_demand = 1.0e-6 * (dissipated_energy + expansion_work) / m_target;  // [MJ/kg]

    if (
        par == nullptr ||
        cpar.target_specie == par->invalid_index ||
        m_target < 10*std::numeric_limits<double>::min()
    )
    {
        energy_demand = SimulationData::infinite_energy_demand;
    }
  
    if (energy_demand < 0.0)
    {
        LOG_ERROR(Error::severity::warning, Error::type::postprocess, "Energy demand is negative: " + std::to_string(energy_demand), cpar.ID);
        energy_demand = SimulationData::infinite_energy_demand;
    }
}


std::string SimulationData::to_csv() const
{
    std::stringstream ss;
    auto format_double = [](std::ostream& os) -> std::ostream& {
        return os << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
    };

    ss << format_double << this->dissipated_energy << ",";
    ss << format_double << this->expansion_work << ",";
    ss << format_double << this->n_target_specie << ",";
    ss << format_double << this->energy_demand << ",";
    ss << format_double << this->R_max << ",";
    ss << format_double << this->T_max << ",";
    ss << format_double << this->t_peak << ",";
    ss << format_double << this->R_min << ",";
    ss << format_double << this->T_min << ",";
    ss << format_double << this->v_max << ",";
    ss << format_double << this->p_internal_max << ",";
    ss << format_double << this->p_internal_min << ",";
    ss << format_double << this->Ma_max << ",";
    ss << format_double << this->T_L_max << ",";
    ss << format_double << this->c_L_max << ",";
    ss << format_double << this->rho_L_max << ",";
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
    ss << std::setw(strw) << "    .expansion_work"     << " = " << format_double << this->expansion_work       << ",    // [J]\n";
    ss << std::setw(strw) << "    .n_target_specie"    << " = " << format_double << this->n_target_specie      << ",    // [mol]\n";
    ss << std::setw(strw) << "    .energy_demand"      << " = " << format_double << this->energy_demand        << ",    // [MJ/kg]\n";
    ss << std::setw(strw) << "    .R_max"              << " = " << format_double << this->R_max                << ",    // [m]\n";
    ss << std::setw(strw) << "    .T_max"              << " = " << format_double << this->T_max                << ",    // [K]\n";
    ss << std::setw(strw) << "    .t_peak"             << " = " << format_double << this->t_peak               << ",    // [s]\n";
    ss << std::setw(strw) << "    .R_min"              << " = " << format_double << this->R_min                << ",    // [m]\n";
    ss << std::setw(strw) << "    .T_min"              << " = " << format_double << this->T_min                << ",    // [K]\n";
    ss << std::setw(strw) << "    .v_max"              << " = " << format_double << this->v_max                << ",    // [m/s]\n";
    ss << std::setw(strw) << "    .p_internal_max"     << " = " << format_double << this->p_internal_max       << ",    // [Pa]\n";
    ss << std::setw(strw) << "    .p_internal_min"     << " = " << format_double << this->p_internal_min       << ",    // [Pa]\n";
    ss << std::setw(strw) << "    .Ma_max"             << " = " << format_double << this->Ma_max               << ",    // [-]\n";
    ss << std::setw(strw) << "    .T_L_max"            << " = " << format_double << this->T_L_max              << ",    // [K]\n";
    ss << std::setw(strw) << "    .c_L_max"            << " = " << format_double << this->c_L_max              << ",    // [m/s]\n";
    ss << std::setw(strw) << "    .rho_L_max"          << " = " << format_double << this->rho_L_max            << ",    // [kg/m^3]\n";
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
    if (ps.ratio->get_num_steps() > 1)
        ss << "ratio=" << format_double << this->cpar.ratio << " [-]; ";
    if (ps.P_amb->get_num_steps() > 1)
        ss << "P_amb=" << format_double << this->cpar.P_amb << " Pa; ";
    if (ps.T_inf->get_num_steps() > 1)
        ss << "T_inf=" << format_double << this->cpar.T_inf << " K; ";
    if (ps.alpha_M->get_num_steps() > 1)
        ss << "alpha_M=" << format_double << this->cpar.alpha_M << "; ";
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
    if (ps.excitation_cycles->get_num_steps() > 1)
        ss << "excitation_cycles=" << format_double << this->cpar.excitation_cycles << "; ";
    if (ps.ramp_up_cycles->get_num_steps() > 1)
        ss << "ramp_up_cycles=" << format_double << this->cpar.ramp_up_cycles << "; ";

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


nlohmann::ordered_json SimulationData::to_json() const
{
    nlohmann::ordered_json j;
    nlohmann::ordered_json postproc;
    postproc["dissipated_energy"] = this->dissipated_energy;
    postproc["expansion_work"] = this->expansion_work;
    postproc["n_target_specie"] = this->n_target_specie;
    postproc["energy_demand"] = this->energy_demand;
    postproc["R_max"] = this->R_max;
    postproc["T_max"] = this->T_max;
    postproc["t_peak"] = this->t_peak;
    postproc["R_min"] = this->R_min;
    postproc["T_min"] = this->T_min;
    postproc["v_max"] = this->v_max;
    postproc["p_internal_max"] = this->p_internal_max;
    postproc["p_internal_min"] = this->p_internal_min;
    postproc["Ma_max"] = this->Ma_max;
    postproc["T_L_max"] = this->T_L_max;
    postproc["c_L_max"] = this->c_L_max;
    postproc["rho_L_max"] = this->rho_L_max;
    j["postproc"] = postproc;
    j["cpar"] = this->cpar.to_json();
    j["sol"] = this->sol.to_json();
    j["excitation"] = nlohmann::ordered_json::object({
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

    const Parameters* par = cpar.par;
    if (par == nullptr) return j;
    j["mechanism"] = nlohmann::ordered_json::object({
        {"model", par->mechanism_name},
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
    nlohmann::ordered_json j = this->to_json();
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

    // Save sol.p_excitation (1D array)
    if (this->sol.p_excitation.size() != this->sol.x.size())
    {
        LOG_ERROR("Mismatch between sol.p_excitation.size() and sol.x.size(): " + std::to_string(this->sol.p_excitation.size()) + " != " + std::to_string(this->sol.x.size()));
        return;
    }
    binary_output_file.write(reinterpret_cast<const char*>(this->sol.p_excitation.data()), this->sol.p_excitation.size() * sizeof(double));

    // Save sol.p_internal (1D array)
    if (this->sol.p_internal.size() != this->sol.x.size())
    {
        LOG_ERROR("Mismatch between sol.p_internal.size() and sol.x.size(): " + std::to_string(this->sol.p_internal.size()) + " != " + std::to_string(this->sol.x.size()));
        return;
    }
    binary_output_file.write(reinterpret_cast<const char*>(this->sol.p_internal.data()), this->sol.p_internal.size() * sizeof(double));

    binary_output_file.close();
}


std::ostream &operator<<(std::ostream &os, const SimulationData &data)
{
    os << data.to_string();
    return os;
}
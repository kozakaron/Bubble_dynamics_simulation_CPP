#include <sstream>
#include <iomanip>

#include "nlohmann/json.hpp"
#include "ode_solver.h"

using json = nlohmann::json;


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


void OdeSolution::push_t_x(const double t_i, const double *x_i)
{
    t.push_back(t_i);
    x.push_back(std::vector<double>(x_i, x_i + num_dim));
}

void OdeSolution::clear()
{
    t.clear();
    x.clear();
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


json OdeSolution::to_json() const
{
    json j;
    j["success"] = this->success();
    j["error"] = ErrorHandler::get_error(this->error_ID).to_string();
    j["runtime"] = this->runtime;
    j["num_steps"] = this->num_steps;
    j["num_saved_steps"] = this->x.size();
    j["num_repeats"] = this->num_repeats;
    j["num_fun_evals"] = this->num_fun_evals;
    j["num_fun_evals_jac"] = this->num_fun_evals_jac;
    j["num_jac_evals"] = this->num_jac_evals;
    j["num_lin_iters"] = this->num_lin_iters;
    j["num_nonlin_iters"] = this->num_nonlin_iters;
    j["total_error"] = this->total_error;
    j["t"] = std::vector<double>({this->t.front(), this->t.back()});
    j["x"] = std::vector<std::vector<double>>({this->x.front(), this->x.back()});
    
    return j;
}



std::ostream &operator<<(std::ostream &os, const OdeSolution &ode)
{
    os << ode.to_string();
    return os;
}
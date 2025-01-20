#include "ode_solver.h"
#include <sstream>
#include <iomanip>


OdeSolution::OdeSolution():
    t(),
    x(),
    num_dim(0),
    num_steps(0),
    num_fun_evals(0),
    num_fun_evals_jac(0),
    num_jac_evals(0),
    num_plu(0),
    num_solve_with_plu(0),
    runtime(0.0),
    error_ID(ErrorHandler::no_error)
{

}


OdeSolution::~OdeSolution()
{ }


is_success OdeSolution::success() const
{
    return error_ID == ErrorHandler::no_error;
}


void OdeSolution::push_t_x(const double t_i, const double *x_i)
{
    t.push_back(t_i);
    x.push_back(std::vector<double>(x_i, x_i + num_dim));
}


std::string percent(const size_t num, const size_t num_steps)
{
    std::stringstream ss;
    ss << " (";
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
    ss << ")";
    return ss.str();
}


std::string OdeSolution::to_string() const
{
    std::stringstream ss;
    size_t strw = 20;
    size_t intw = 10;
    
    ss << std::setw(strw) << "success: ";
    if (success())
    {
        ss << colors::green << colors::bold << std::setw(intw) << "true" << colors::reset << "\n";
    }
    else
    {
        ss << colors::red << colors::bold << std::setw(intw) << "false" << colors::reset << "\n";
        ss << std::setw(strw) << "error: " << ErrorHandler::get_error(error_ID) << "\n";
    }

    ss << std::setw(strw) << "runtime: "            << std::setw(intw+3) << Timer::format_time(runtime) << "\n";
    ss << std::setw(strw) << "num_steps: "          << std::setw(intw) << num_steps          << percent(num_steps, num_steps) << "\n";
    ss << std::setw(strw) << "num_fun_evals: "      << std::setw(intw) << num_fun_evals      << percent(num_fun_evals, num_steps) << "\n";
    ss << std::setw(strw) << "num_fun_evals_jac: "  << std::setw(intw) << num_fun_evals_jac  << percent(num_fun_evals_jac, num_steps) << "\n";
    ss << std::setw(strw) << "num_jac_evals: "      << std::setw(intw) << num_jac_evals      << percent(num_jac_evals, num_steps) << "\n";
    ss << std::setw(strw) << "num_plu: "            << std::setw(intw) << num_plu            << percent(num_plu, num_steps) << "\n";
    ss << std::setw(strw) << "num_solve_with_plu: " << std::setw(intw) << num_solve_with_plu << percent(num_solve_with_plu, num_steps) << "\n";
    
    if (this->t.size() >=2 && this->x.size() >= 2 && this->num_dim > 0)
    {
        ss << std::setw(strw) << "t: " << "{";
        ss << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10) << this->t[0]                  << ", ..., ";
        ss << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10) << this->t[this->t.size() - 1] << "}\n";
        ss << std::setw(strw) << "x: " << "{\n";
        ss << std::setw(strw+4) << " " << ::to_string((double*)this->x[0].data(),                  this->x[0].size())                  << ",\n";
        ss << std::setw(strw+4) << " " << " ..., \n";
        ss << std::setw(strw+4) << " " << ::to_string((double*)this->x[this->x.size() - 1].data(), this->x[this->x.size() - 1].size()) << "\n";
        ss << std::setw(strw) << " " << "}\n";
    }

    return ss.str();
}


std::ostream &operator<<(std::ostream &os, const OdeSolution &ode)
{
    os << ode.to_string();
    return os;
}
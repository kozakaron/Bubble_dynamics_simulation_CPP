#include "ode_solver.h"
#include <sstream>
#include <iomanip>
#include <cmath>


OdeSolution::OdeSolution():
    t(),
    x(),
    num_dim(0),
    num_steps(0),
    num_repeats(0),
    num_fun_evals(0),
    num_fun_evals_jac(0),
    num_jac_evals(0),
    num_plu(0),
    num_solve_with_plu(0),
    total_error(0.0),
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
    num_plu = 0;
    num_solve_with_plu = 0;
    total_error = 0.0;
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
    ss << this->num_plu << "," << this->num_solve_with_plu << "," << format_double << this->total_error << ",";
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
    ss << std::setw(strw) << "    .total_error"        << " = " << std::setw(intw) << std::scientific << total_error << ",\n";
    ss << std::setw(strw) << "    .num_steps"          << " = " << std::setw(intw) << std::to_string(num_steps)          + "," << percent(num_steps, num_steps) << "\n";
    ss << std::setw(strw) << "    .num_saved_steps"    << " = " << std::setw(intw) << std::to_string(x.size())           + "," << percent(x.size(), num_steps) << "\n";
    ss << std::setw(strw) << "    .num_repeats"        << " = " << std::setw(intw) << std::to_string(num_repeats)        + "," << percent(num_repeats, num_steps) << "\n";
    ss << std::setw(strw) << "    .num_fun_evals"      << " = " << std::setw(intw) << std::to_string(num_fun_evals)      + "," << percent(num_fun_evals, num_steps) << "\n";
    ss << std::setw(strw) << "    .num_fun_evals_jac"  << " = " << std::setw(intw) << std::to_string(num_fun_evals_jac)  + "," << percent(num_fun_evals_jac, num_steps) << "\n";
    ss << std::setw(strw) << "    .num_jac_evals"      << " = " << std::setw(intw) << std::to_string(num_jac_evals)      + "," << percent(num_jac_evals, num_steps) << "\n";
    ss << std::setw(strw) << "    .num_plu"            << " = " << std::setw(intw) << std::to_string(num_plu)            + "," << percent(num_plu, num_steps) << "\n";
    ss << std::setw(strw) << "    .num_solve_with_plu" << " = " << std::setw(intw) << std::to_string(num_solve_with_plu) + "," << percent(num_solve_with_plu, num_steps) << "\n";
    
    if (this->t.size() >=2 && this->x.size() >= 2 && this->num_dim > 0)
    {
        ss << std::setw(strw) << "    .t" << " = " << "{";
        ss << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10) << this->t[0]                  << ", /* ..., */ ";
        ss << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10) << this->t[this->t.size() - 1] << "},\n";
        ss << std::setw(strw) << "    .x" << " = " << "{\n";
        ss << std::setw(strw+8) << " " << ::to_string((double*)this->x[0].data(),                  this->x[0].size())                  << ",\n";
        ss << std::setw(strw+8) << " " << " // ..., \n";
        ss << std::setw(strw+8) << " " << ::to_string((double*)this->x[this->x.size() - 1].data(), this->x[this->x.size() - 1].size()) << "\n";
        ss << std::setw(strw+4) << " " << "}\n";
    }

    if (with_code) ss << "}";
    ss << std::right;
    return ss.str();
}


std::ostream &operator<<(std::ostream &os, const OdeSolution &ode)
{
    os << ode.to_string();
    return os;
}


OdeSolver::OdeSolver(
    const double rel_tol,
    const double abs_tol,
    const double dt_min,
    const double dt_max
):
    error_ID(nullptr),
    x(nullptr),
    x_new(nullptr),
    loc_error(nullptr),
    rel_tol(rel_tol),
    abs_tol(abs_tol),
    dt_min(dt_min),
    dt_max(dt_max)
{ }


OdeSolver::~OdeSolver()
{
    if (this->x != nullptr) delete[] this->x;
    if (this->x_new != nullptr) delete[] this->x_new;
    if (this->loc_error != nullptr) delete[] this->loc_error;
}


double norm(const double *x, const size_t num_dim)
{
    double norm = 0.0;
    for (size_t i = 0; i < num_dim; ++i)
    {
        norm += x[i] * x[i];
    }
    return std::sqrt(norm);
}


is_success OdeSolver::preprocess(
    const double t_int_0,
    const double t_int_1,
    const double *x_0,
    const size_t num_dim,
    const OdeFun_t& ode_fun,
    size_t *error_ID
) {
    this->timer.start();
    this->sol.clear();
    if (t_int_1 <= t_int_0)
    {
        std::string message = "Time interval is invalid: t_int_1 <= t_int_0 (" + std::to_string(t_int_1) + " <= " + std::to_string(t_int_0) + ")";
        this->sol.error_ID = LOG_ERROR(Error::severity::error, Error::type::preprocess, message);
        return false;
    } else {
        this->t_int = std::make_pair(t_int_0, t_int_1);
    }
    bool change_num_dim = this->sol.num_dim != num_dim;
    this->sol.num_dim = num_dim;
    this->ode_fun = ode_fun;
    this->error_ID = error_ID;

    if (change_num_dim || this->x == nullptr || this->x_new == nullptr || this->loc_error == nullptr)
    {
        if (this->x != nullptr) delete[] this->x;
        if (this->x_new != nullptr) delete[] this->x_new;
        if (this->loc_error != nullptr) delete[] this->loc_error;
        this->x = new double[num_dim];
        this->x_new = new double[num_dim];
        this->loc_error = new double[num_dim];
    }
    if (x_0 != nullptr)
        std::copy(x_0, x_0 + num_dim, this->x);
    this->sol.push_t_x(t_int_0, x_0);

    this->t = t_int_0;
    /*if (!this->ode_fun(this->t_int.first, this->x, this->x_new))
    {
        return false;
    }
    this->dt = 0.01 * (this->abs_tol + this->rel_tol * norm(this->x, num_dim)) /\
                      (this->abs_tol + this->rel_tol * norm(this->x_new, num_dim));
    this->dt = std::max(this->dt, this->dt_min);
    this->dt = std::min(this->dt, this->dt_max);
    this->dt = std::min(this->dt, 0.01 * (this->t_int.second - this->t_int.first));*/
    this->dt = 1.0e-20;

    return true;
}


is_success OdeSolver::postprocess(bool save_solution)
{
    if (this->error_ID != nullptr && this->sol.error_ID == ErrorHandler::no_error)
    {
        this->sol.error_ID = *this->error_ID;
    }
    if (!save_solution)
    {
        this->sol.push_t_x(this->t, this->x);
    }
    this->sol.runtime = this->timer.lap();

    return this->sol.success();
}


double OdeSolver::step_size_control(
    const double loc_error,
    const double tol
) const {
    if (!std::isfinite(loc_error))
        return this->dt_min;
    if (this->t >= this->t_int.second)
        return this->dt;

    const double loc_error_nonzero = std::max(loc_error, 1.0e-30);
    double factor = 0.9 * std::pow(tol / loc_error_nonzero, 1.0 / (this->order + 1));
    //factor = std::max(0.1, std::min(10.0, factor));
    double dt_new = this->dt * factor;
    dt_new = std::max(this->dt_min, dt_new);
    dt_new = std::min(this->dt_max, dt_new);
    dt_new = std::min(dt_new, this->t_int.second - this->t);

    if (!std::isfinite(dt_new))
        return this->dt_min;
    return dt_new;
}


is_success OdeSolver::solve(
    const double t_int_0,
    const double t_int_1,
    const double *x_0,
    const size_t num_dim,
    const OdeFun_t& ode_fun,
    size_t *error_ID,
    double timeout,
    bool save_solution
) {
    if (!this->preprocess(t_int_0, t_int_1, x_0, num_dim, ode_fun, error_ID))
    {
        return false;
    }

    bool success = true;
    while (this->t < this->t_int.second)
    {
    // Break conditions
        if (!this->step())
        {
            success = false;
            break;
        }
        if (this->error_ID != nullptr && *this->error_ID != ErrorHandler::no_error)
        {
            success = false;
            break;
        }
        if (this->timer.lap() > timeout)
        {
            this->sol.error_ID = LOG_ERROR(Error::severity::error, Error::type::timeout, "Timeout: " + std::to_string(timeout) + " s", 0);
            success = false;
            break;
        }
    // Accept/reject step
        const double loc_error = norm(this->loc_error, this->sol.num_dim);
        const double tol = this->abs_tol + this->rel_tol * norm(this->x_new, this->sol.num_dim);
        if ((success && loc_error < tol) || this->dt <= this->dt_min)
        {   
            // Accept step
            this->t = this->t + this->dt;
            std::swap(this->x, this->x_new);
            this->sol.num_steps++;
            this->sol.total_error += loc_error;
            if (save_solution)
            {
                this->sol.push_t_x(this->t, this->x);
            }
        } else {
            // Reject step
            this->sol.num_repeats += 1;
        }
        
        this->dt = this->step_size_control(loc_error, tol);
    }

    return this->postprocess(save_solution);
}


OdeSolution& OdeSolver::get_solution() const
{
    return const_cast<OdeSolution&>(this->sol);
}


RKCK45::RKCK45(
    const double rel_tol,
    const double abs_tol,
    const double dt_min,
    const double dt_max
):
    OdeSolver(rel_tol, abs_tol, dt_min, dt_max),
    k(new double*[6]),
    x_stage(nullptr)
{
    this->order = 4;
    this->stages = 6;
    for (size_t i = 0; i < this->stages; ++i)
    {
        this->k[i] = nullptr;
    }
}


RKCK45::~RKCK45()
{
    for (size_t i = 0; i < this->stages; ++i)
    {
        if (this->k[i] != nullptr) delete[] this->k[i];
    }
    delete[] this->k;

    if (this->x_stage != nullptr) delete[] this->x_stage;
}


is_success RKCK45::preprocess(
    const double t_int_0,
    const double t_int_1,
    const double* x_0,
    const size_t num_dim,
    const OdeFun_t& ode_fun,
    size_t* error_ID
) {
    is_success success = OdeSolver::preprocess(t_int_0, t_int_1, x_0, num_dim, ode_fun, error_ID);

    bool change_num_dim = this->sol.num_dim != num_dim;
    if (change_num_dim || this->k[0] == nullptr || this->x_stage == nullptr)
    {
        for (size_t i = 0; i < this->stages; ++i)
        {
            if (this->k[i] != nullptr) delete[] this->k[i];
            this->k[i] = new double[num_dim];
        }
        if (this->x_stage != nullptr) delete[] this->x_stage;
        this->x_stage = new double[num_dim];
    }

    return success;
}


is_success RKCK45::step()
{
    is_success success = true;
    for (size_t stage = 0; stage < this->stages; ++stage)
    {
        const double t_stage = this->t + this->dt * RKCK45::c[stage];
        for (size_t i = 0; i < this->sol.num_dim; ++i)
        {
            this->x_stage[i] = 0.0;
            for (size_t j = 0; j < stage; ++j)
            {
                this->x_stage[i] += RKCK45::a[stage][j] * this->k[j][i];
            }
            this->x_stage[i] = this->x[i] + this->dt * this->x_stage[i];
        }
        success &= this->ode_fun(t_stage, this->x_stage, this->k[stage]);
        if (!success)
        {
            break;
        }
    }

    if (success) {
        for (size_t i = 0; i < this->sol.num_dim; ++i)
        {
            this->x_new[i] = 0.0;
            this->loc_error[i] = 0.0;
            for (size_t stage = 0; stage < this->stages; ++stage)
            {
                this->x_new[i] += RKCK45::b_lower[stage] * this->k[stage][i];
                this->loc_error[i] += RKCK45::b_error[stage] * this->k[stage][i];
            }
            this->x_new[i] = this->x[i] + this->dt * this->x_new[i];
            this->loc_error[i] = this->dt * this->loc_error[i];
        }
    }

    return success;
}
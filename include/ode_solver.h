#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H
#include <vector>
#include <ostream>

#include "nlohmann/json_fwd.hpp"
#include "common.h"
#include "ode_fun.h"

class OdeSolution
{
public:
    static constexpr char csv_header[] = "success,num_dim,num_steps,num_repeats,num_fun_evals,num_fun_evals_jac,num_jac_evals,num_lin_iters,num_nonlin_iters,runtime,t_last,R_max,T_max,t_max,x_0,x_last";

    std::vector<double> t;
    std::vector<std::vector<double>> x;
    double R_max;
    double T_max;
    double t_max;
    size_t num_dim;
    size_t num_steps;
    size_t num_repeats;
    size_t num_fun_evals;
    size_t num_fun_evals_jac;
    size_t num_jac_evals;
    size_t num_lin_iters;
    size_t num_nonlin_iters;
    std::vector<double> total_error;
    
    double runtime;
    size_t error_ID;


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


class OdeSolver
{
public:
    virtual ~OdeSolver() = default;
    virtual OdeSolution solve(
        const double t_max,
        OdeFun* ode_ptr,
        double timeout = 1.0e30,
        bool save_solution = false,
        bool save_jacobian = false
    ) = 0;
};



#endif // ODE_SOLVER_H
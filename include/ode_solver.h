#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H
#include <vector>
#include <ostream>

#include "common.h"
#include "ode_fun.h"

class OdeSolution
{
public:
    static constexpr char csv_header[] = "success,num_dim,num_steps,num_repeats,num_fun_evals,num_fun_evals_jac,num_jac_evals,num_plu,num_solve_with_plu,total_error,runtime,t_last,x_0,x_last";

    std::vector<double> t;
    std::vector<std::vector<double>> x;
    size_t num_dim;
    size_t num_steps;
    size_t num_repeats;
    size_t num_fun_evals;
    size_t num_fun_evals_jac;
    size_t num_jac_evals;
    size_t num_plu; // TODO: change
    size_t num_solve_with_plu;
    double total_error;
    
    double runtime;
    size_t error_ID;


    OdeSolution();
    ~OdeSolution();
    is_success success() const;
    void push_t_x(const double t_i, const double *x_i);
    void clear();
    std::string to_csv() const;
    std::string to_string(const bool colored=true, const bool with_code=true) const;
    friend std::ostream &operator<<(std::ostream &os, const OdeSolution &ode);
};


class OdeSolver
{
public:
    virtual OdeSolution solve(
        const double t_max,
        OdeFun* ode_ptr,
        double timeout = 1.0e30,
        bool save_solution = false
    ) = 0;
};



#endif // ODE_SOLVER_H
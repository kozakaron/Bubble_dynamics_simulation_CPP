#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H
#include <vector>
#include <ostream>

#include "common.h"

class OdeSolution
{
public:
    std::vector<double> t;
    std::vector<std::vector<double>> x;
    size_t num_dim;
    size_t num_steps;
    size_t num_fun_evals;
    size_t num_fun_evals_jac;
    size_t num_jac_evals;
    size_t num_plu;
    size_t num_solve_with_plu;
    
    double runtime;
    size_t error_ID;


    OdeSolution();
    ~OdeSolution();
    is_success success() const;
    void push_t_x(const double t_i, const double *x_i);
    std::string to_string() const;
    friend std::ostream &operator<<(std::ostream &os, const OdeSolution &ode);
};



#endif // ODE_SOLVER_H
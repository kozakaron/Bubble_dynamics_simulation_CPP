#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H
#include <vector>
#include <ostream>
#include <functional>

#include "common.h"

class OdeSolution
{
public:
    std::vector<double> t;
    std::vector<std::vector<double>> x;
    size_t num_dim;
    size_t num_steps;
    size_t num_repeats;
    size_t num_fun_evals;
    size_t num_fun_evals_jac;
    size_t num_jac_evals;
    size_t num_plu;
    size_t num_solve_with_plu;
    double total_error;
    
    double runtime;
    size_t error_ID;


    OdeSolution();
    ~OdeSolution();
    is_success success() const;
    void push_t_x(const double t_i, const double *x_i);
    void clear();
    std::string to_string() const;
    friend std::ostream &operator<<(std::ostream &os, const OdeSolution &ode);
};


class OdeSolver
{
public:
    using OdeFun_t = std::function<is_success(const double, const double*, double*)>;

    OdeSolution sol;
    Timer timer;
    OdeFun_t ode_fun;
    size_t* error_ID;
    size_t order;   // TODO: seperate order of embedded method
    size_t stages;

    std::pair<double, double> t_int;
    double *x;
    double *x_new;
    double *loc_error;
    double t;
    const double rel_tol;
    const double abs_tol;
    const double dt_min;
    const double dt_max;
    double dt;

    OdeSolver(
        const double rel_tol=1.0e-10,
        const double abs_tol=1.0e-10,
        const double dt_min=1.0e-20,
        const double dt_max=1.0e-1
    );
    ~OdeSolver();
    virtual is_success preprocess(
        const double t_int_0,
        const double t_int_1,
        const double *x_0,
        const size_t num_dim,
        const OdeFun_t& ode_fun,
        size_t *error_ID
    );
    is_success postprocess(
        bool save_solution
    );
    virtual double step_size_control(
        const double loc_error,
        const double tol
    ) const;
    virtual is_success step() = 0;
    is_success solve(
        const double t_int_0,
        const double t_int_1,
        const double *x_0,
        const size_t num_dim,
        const OdeFun_t& ode_fun,
        size_t *error_ID,
        double timeout=1.0e30,
        bool save_solution=false
    );
    OdeSolution& get_solution() const;
};


class RKCK45 : public OdeSolver
{
public:
    static constexpr double a[6][6] = {
        {0.0,            0.0,          0.0,            0.0,               0.0},
        {1.0/5.0,        0.0,          0.0,            0.0,               0.0},
        {3.0/40.0,       9.0/40.0,     0.0,            0.0,               0.0},
        {3.0/10.0,      -9.0/10.0,     6.0/5.0,        0.0,               0.0},
        {-11.0/54.0,     5.0/2.0,     -70.0/27.0,      35.0/27.0,         0.0},
        {1631.0/55296.0, 175.0/512.0,  575.0/13824.0,  44275.0/110592.0,  253.0/4096.0}
     };
    static constexpr double b_error[] = {-0.0042937748015873,  0.                ,  0.0186685860938579, -0.0341550268308081, -0.0193219866071429,  0.0391022021456804};
    static constexpr double b_lower[] = {2825.0/27648.0, 0.0, 18575.0/48384.0, 13525.0/55296.0, 277.0/14336.0, 1.0/4.0};
    static constexpr double c[] = {0.0, 1.0/5.0, 3.0/10.0, 3.0/5.0, 1.0, 7.0/8.0};

    double **k;
    double *x_stage;

    RKCK45(
        const double rel_tol=1.0e-10,
        const double abs_tol=1.0e-10,
        const double dt_min=1.0e-30,
        const double dt_max=1.0e30
    );
    ~RKCK45();
    is_success preprocess(
        const double t_int_0,
        const double t_int_1,
        const double* x_0,
        const size_t num_dim,
        const OdeFun_t& ode_fun,
        size_t* error_ID
    ) override;
    is_success step() override;
};

#endif // ODE_SOLVER_H
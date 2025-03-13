#ifndef ODE_SOLVER_SUNDIALS_H
#define ODE_SOLVER_SUNDIALS_H

#include <cvode/cvode.h>                // prototypes for CVODE fcts., consts.
#include <nvector/nvector_serial.h>     // access to serial N_Vector
#include <sunlinsol/sunlinsol_dense.h>  // access to dense SUNLinearSolver
#include <sunmatrix/sunmatrix_dense.h>  // access to dense SUNMatrix

#include "common.h"
#include "ode_fun.h"
#include "ode_solver.h"

class OdeSolverCVODE
{
public:
    size_t init_error_ID;   // only used to indicate errors in the constructor
    SUNContext sun_context;
    double t;
    N_Vector x;
    N_Vector abstol;
    SUNMatrix A;
    SUNLinearSolver linear_solver;
    void* cvode_mem;

    OdeSolverCVODE(const size_t num_dim);
    ~OdeSolverCVODE();
    OdeSolution solve(
        const double t_max,
        OdeFun* ode,
        double timeout = 1.0e30,
        bool save_solution = false
    );
};


#endif // ODE_SOLVER_SUNDIALS_H
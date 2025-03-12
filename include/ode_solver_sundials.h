#ifndef ODE_SOLVER_SUNDIALS_H
#define ODE_SOLVER_SUNDIALS_H

#include <cvode/cvode.h>                // prototypes for CVODE fcts., consts.
#include <nvector/nvector_serial.h>     // access to serial N_Vector
#include <sunlinsol/sunlinsol_dense.h>  // access to dense SUNLinearSolver
#include <sunmatrix/sunmatrix_dense.h>  // access to dense SUNMatrix

#include "common.h"
#include "ode_fun.h"

class OdeSolverCVODE
{
public:
    size_t error_ID;
    SUNContext sun_context;
    double t;
    N_Vector x;
    N_Vector abstol;
    SUNMatrix A;
    SUNLinearSolver linear_solver;
    void* cvode_mem;

    OdeSolverCVODE(OdeFun *ode);
    ~OdeSolverCVODE();
};


#endif // ODE_SOLVER_SUNDIALS_H
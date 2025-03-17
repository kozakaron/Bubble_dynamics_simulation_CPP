#ifndef ODE_SOLVER_SUNDIALS_H
#define ODE_SOLVER_SUNDIALS_H

#include <cvode/cvode.h>                // prototypes for CVODE fcts., consts.
#include <nvector/nvector_serial.h>     // access to serial N_Vector
#include <sunlinsol/sunlinsol_dense.h>  // access to dense SUNLinearSolver
#include <sunmatrix/sunmatrix_dense.h>  // access to dense SUNMatrix

#include "common.h"
#include "ode_fun.h"
#include "ode_solver.h"

// Used to pass data to the (CVRhsFn) right_hand_side function as void* user_data. 
// Used to pass data to the (SUNErrHandlerFn) error_function function as void* err_user_data.
struct UserData
{
    // right_hand_side function arguments
    OdeFun* ode_ptr      = nullptr;
    Timer* timer_ptr     = nullptr;
    double timeout       = 1.0e30;
    bool timed_out       = false;

    // error_function function arguments
    size_t* error_ID_ptr = nullptr;
    size_t ID            = 0;
};


class OdeSolverCVODE
{
public:
    SUNContext sun_context;           // An opaque pointer used by SUNDIALS objects for error handling, logging, profiling, etc.
    double t;                         // simulation time
    N_Vector x;                       // simulation state vector
    N_Vector abstol;                  // absolute tolerance vector
    SUNMatrix A;                      // matrix for linear solver
    SUNLinearSolver linear_solver;    // linear solver
    void* cvode_mem;                  // CVODE memory block
    size_t init_error_ID;             // error ID from the initialization (solution.error_ID is used during simulation)
    UserData user_data;               // to pass data to the right_hand_side function and error_function

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
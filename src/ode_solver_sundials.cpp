#include "ode_solver_sundials.h"

#include <type_traits>
static_assert(std::is_same<sunrealtype, double>::value, "sunrealtype must be double");

// Macro for handling SUNDIALS return codes.
// Use inside a function as: HANDLE_ERROR_CODE(CVodeCreate(...));
// instead of: int retval = CVodeCreate(...); if (retval != 0) { ... }
#define HANDLE_ERROR_CODE(...) \
{ \
    int retval = __VA_ARGS__; \
    if (retval != 0) \
    { \
        *error_ID = LOG_ERROR(Error::severity::error, Error::type::cvode, "CVODE init error: " # __VA_ARGS__ "returned with code " + std::to_string(retval)); \
        return; \
    } \
}

// Macro for handling pointers returned by SUNDIALS.
// Use inside a function as: HANDLE_RETURN_PTR(pointer, N_VNew_Serial(...));
// instead of: pointer = N_VNew_Serial(...); if (pointer == nullptr) { ... }
#define HANDLE_RETURN_PTR(pointer, ...) \
{ \
    pointer = __VA_ARGS__; \
    if (pointer == nullptr) \
    { \
        *error_ID = LOG_ERROR(Error::severity::error, Error::type::cvode, "CVODE init error: " # __VA_ARGS__ "returned with nullptr"); \
        return; \
    } \
}


// CVRhsFn, CVMonitorFn, SUNErrHandlerFn definitisions

struct UserData
{
    OdeFun* ode;
    Timer* timer;
    const double timeout;
    bool timed_out = false;
};


struct ErrorUserData
{
    size_t* error_ID;
    size_t ID = 0;  // TODO: save cpar ID in errors
};


int right_hand_side(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
    UserData* data = (UserData*)user_data;
    OdeFun* ode = data->ode;
    Timer* timer = data->timer;
    if (timer->lap() > data->timeout)
    {
        data->timed_out = true;
        return 1;
    }
    is_success ret = ode->operator()(t, NV_DATA_S(y), NV_DATA_S(ydot));
    return !ret;
}
static_assert(std::is_same<decltype(&right_hand_side), CVRhsFn>::value, "right_hand_side must match CVRhsFn");


void error_function(int line, const char* func, const char* file, const char* msg, SUNErrCode err_code, void* err_user_data, SUNContext sunctx)
{
    size_t *error_ID = (size_t*)err_user_data;
    Error error(
        Error::severity::error,
        Error::type::cvode,
        "CVODE error code " + std::to_string(err_code) + ": " + std::string(msg),
        std::string(func),
        std::string(file),
        line
    );
    *error_ID = ErrorHandler::log_error(error);
}
static_assert(std::is_same<decltype(&error_function), SUNErrHandlerFn>::value, "error_function must match SUNErrHandlerFn");


// OdeSolverCVODE definitions

OdeSolverCVODE::OdeSolverCVODE(const size_t num_dim):
    sun_context(nullptr),
    t(0.0),
    x(nullptr),
    abstol(nullptr),
    A(nullptr),
    linear_solver(nullptr),
    cvode_mem(nullptr),
    init_error_ID(ErrorHandler::no_error)
{
    const double _abstol = 1e-10;
    const double _reltol = 1e-10;

    // Create the SUNDIALS context
    size_t* error_ID = &(this->init_error_ID);  // errors are saved in this->init_error_ID
    HANDLE_ERROR_CODE(SUNContext_Create(SUN_COMM_NULL, &sun_context));
    HANDLE_ERROR_CODE(SUNContext_ClearErrHandlers(sun_context));
    HANDLE_ERROR_CODE(SUNContext_PushErrHandler(sun_context, error_function, error_ID));

    // Setup vectors
    HANDLE_RETURN_PTR(x, N_VNew_Serial(num_dim, sun_context));
    /*HANDLE_RETURN_PTR(abstol, N_VNew_Serial(num_dim, sun_context))
    for (int i = 0; i < num_dim; i++)
        NV_Ith_S(abstol, i) = _abstol;
    HANDLE_ERROR_CODE(CVodeSVtolerances(cvode_mem, _reltol, abstol));
    */ // TODO: set abstol elementwise for better performance

    // Setup CVODE
    HANDLE_RETURN_PTR(cvode_mem, CVodeCreate(CV_BDF, sun_context));
    HANDLE_ERROR_CODE(CVodeInit(cvode_mem, right_hand_side, 0.0, x));
    HANDLE_ERROR_CODE(CVodeSetMaxNumSteps(cvode_mem, 10000000));
    HANDLE_ERROR_CODE(CVodeSStolerances(cvode_mem, _reltol, _abstol));

    // Setup linear solver
    HANDLE_RETURN_PTR(A, SUNDenseMatrix(num_dim, num_dim, sun_context));
    HANDLE_RETURN_PTR(linear_solver, SUNLinSol_Dense(x, A, sun_context));
    // TODO: investigate nonlinear solvers
    HANDLE_ERROR_CODE(CVodeSetLinearSolver(cvode_mem, linear_solver, A));
}


OdeSolverCVODE::~OdeSolverCVODE()
{
    N_VDestroy(x);
    N_VDestroy(abstol);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(linear_solver);
    SUNMatDestroy(A);
    SUNContext_Free(&sun_context);
}


// used in OdeSolverCVODE::solve
void init_solve(
    SUNContext sun_context,
    void* cvode_mem,
    void* user_data,
    N_Vector x,
    size_t* error_ID
)
{
    // errors are saved in solution.error_ID
    HANDLE_ERROR_CODE(SUNContext_ClearErrHandlers(sun_context));
    HANDLE_ERROR_CODE(SUNContext_PushErrHandler(sun_context, error_function, error_ID));
    HANDLE_ERROR_CODE(CVodeReInit(cvode_mem, 0.0, x));
    HANDLE_ERROR_CODE(CVodeSetUserData(cvode_mem, user_data));
}


// used in OdeSolverCVODE::solve
void construct_solution(
    OdeSolution& solution,
    void* cvode_mem,
    size_t* error_ID
)
{
    long int helper1, helper2;
    HANDLE_ERROR_CODE(CVodeGetNumSteps(cvode_mem, &helper1));
    solution.num_steps = helper1;
    HANDLE_ERROR_CODE(CVodeGetNumRhsEvals(cvode_mem, &helper1));
    solution.num_fun_evals = helper1;
    HANDLE_ERROR_CODE(CVodeGetNumErrTestFails(cvode_mem, &helper1));
    HANDLE_ERROR_CODE(CVodeGetNumStepSolveFails(cvode_mem, &helper2));
    solution.num_repeats = helper1 + helper2;
    HANDLE_ERROR_CODE(CVodeGetNumJacEvals(cvode_mem, &helper1));
    HANDLE_ERROR_CODE(CVodeGetNumJtimesEvals(cvode_mem, &helper2));
    solution.num_jac_evals = helper1 + helper2;
    HANDLE_ERROR_CODE(CVodeGetNumLinRhsEvals(cvode_mem, &helper1));
    solution.num_fun_evals_jac = helper1;
}


OdeSolution OdeSolverCVODE::solve(
    const double t_max,
    OdeFun* ode,
    double timeout,
    bool save_solution
)
{
    // Setup
    Timer timer; timer.start();
    OdeSolution solution;
    auto user_data = UserData{
        .ode     = ode,
        .timer   = &timer,
        .timeout = timeout
    };

    if (this->init_error_ID != ErrorHandler::no_error)
    {
        solution.error_ID = this->init_error_ID;
        return solution;
    }

    solution.num_dim = NV_LENGTH_S(x);
    ode->initial_conditions(NV_DATA_S(x));
    solution.push_t_x(0.0, NV_DATA_S(x));
    init_solve(sun_context, cvode_mem, &user_data, x, &(solution.error_ID));
    if (solution.error_ID != ErrorHandler::no_error)  return solution;

    // Solve
    int itask = save_solution ? CV_ONE_STEP : CV_NORMAL;
    while (true)
    {
        // Integration
        int retval = CVode(cvode_mem, t_max, x, &t, itask);

        // Success
        if (retval == CV_SUCCESS)
        {
            solution.push_t_x(t, NV_DATA_S(x));
        }

        // Failure
        else
        {
            // TODO: Build this into error_function, mention simulation time t
            if (user_data.timed_out)
            {
                solution.error_ID = LOG_ERROR(
                    Error::severity::error, Error::type::timeout,
                    "CVODE error: CVode timed out after " + std::to_string(timeout) + " s",
                    ode->cpar.ID
                );
            }
        }

        // Exit conditions
        if (t >= t_max)  break;
        if (user_data.timed_out)  break;
        if (retval != CV_SUCCESS)  break;
    }

    // fill solution
    construct_solution(solution, cvode_mem, &(solution.error_ID));
    solution.runtime = timer.lap();

    return solution;
}
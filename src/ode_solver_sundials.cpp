#include "ode_solver_sundials.h"

#include <type_traits>
static_assert(std::is_same<sunrealtype, double>::value, "sunrealtype must be double");

// Macro for handling SUNDIALS return codes.
// Use as: HANDLE_ERROR_CODE(CVodeCreate(...));
// instead of: int retval = CVodeCreate(...); if (retval != 0) { ... }
#define HANDLE_ERROR_CODE(...) \
{ \
    (void*) error_ID_ptr; /* a variable named 'size_t* error_ID_ptr;' must be in the scope of the macro call*/ \
    static_assert(std::is_same<decltype(error_ID_ptr), size_t*>::value, "error_ID_ptr must be size_t*"); \
    \
    int retval = __VA_ARGS__; \
    if (retval != 0) \
    { \
        *error_ID_ptr = LOG_ERROR(Error::severity::error, Error::type::cvode, "CVODE init error: " # __VA_ARGS__ "returned with code " + std::to_string(retval)); \
        return; \
    } \
}

// Macro for handling pointers returned by SUNDIALS.
// Use as: HANDLE_RETURN_PTR(pointer, N_VNew_Serial(...));
// instead of: pointer = N_VNew_Serial(...); if (pointer == nullptr) { ... }
#define HANDLE_RETURN_PTR(pointer, ...) \
{ \
    (void*) error_ID_ptr; /* a variable named 'size_t* error_ID_ptr;' must be in the scope of the macro call*/ \
    static_assert(std::is_same<decltype(error_ID_ptr), size_t*>::value, "error_ID_ptr must be size_t*"); \
    \
    pointer = __VA_ARGS__; \
    if (pointer == nullptr) \
    { \
        *error_ID_ptr = LOG_ERROR(Error::severity::error, Error::type::cvode, "CVODE init error: " # __VA_ARGS__ "returned with nullptr"); \
        return; \
    } \
}


// CVRhsFn, SUNErrHandlerFn definitisions:

// Check if user_data is valid (holds valid pointers).
bool check_user_data(void* user_data)
{
    UserData* data = (UserData*)user_data;
    if (data == nullptr)
    {
        LOG_ERROR("user_data is nullptr, unexpected behaviour");
        return false;
    }

    OdeFun* ode_ptr = data->ode_ptr;
    if (ode_ptr == nullptr)
    {
        LOG_ERROR("ode_ptr is nullptr, unexpected behaviour");
        return false;
    }

    Timer* timer_ptr = data->timer_ptr;
    if (timer_ptr == nullptr)
    {
        LOG_ERROR("timer_ptr is nullptr, unexpected behaviour");
        return false;
    }

    size_t* error_ID_ptr = data->error_ID_ptr;
    if (error_ID_ptr == nullptr)
    {
        LOG_ERROR("error_ID_ptr is nullptr, unexpected behaviour");
        return false;
    }

    return true;
}


// Right hand side function for CVODE (CVRhsFn)
int right_hand_side(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
    // gather user data
    if (!check_user_data(user_data))  return 1;
    UserData* data = (UserData*)user_data;
    OdeFun* ode_ptr = data->ode_ptr;
    Timer* timer_ptr = data->timer_ptr;

    // check timeout
    if (timer_ptr->lap() > data->timeout)
    {
        data->timed_out = true;
        return 1;
    }

    // calculate right hand side
    is_success ret = ode_ptr->operator()(t, NV_DATA_S(y), NV_DATA_S(ydot));    // ydot = f(t, y)
    return !ret;
}
static_assert(std::is_same<decltype(&right_hand_side), CVRhsFn>::value, "right_hand_side must match CVRhsFn");


// Error handling function for CVODE (SUNErrHandlerFn)
void error_function(int line, const char* func, const char* file, const char* msg, SUNErrCode err_code, void* err_user_data, SUNContext sunctx)
{
    (void) sunctx;
    if (!check_user_data(err_user_data))  return;
    UserData* user_data = (UserData*)err_user_data;
    size_t* error_ID_ptr = user_data->error_ID_ptr;
    Error error;
    
    if (user_data->timed_out)
    {
        error = Error(
            Error::severity::error,
            Error::type::timeout,
            "CVODE error: CVode timed out after " + std::to_string(user_data->timeout) + " s",
            std::string(func),
            std::string(file),
            line,
            user_data->ID
        );
    } else {
        error = Error(
            Error::severity::error,
            Error::type::cvode,
            "CVODE error code " + std::to_string(err_code) + ": " + std::string(msg),
            std::string(func),
            std::string(file),
            line,
            user_data->ID
        );
    }

    *error_ID_ptr = ErrorHandler::log_error(error);
}
static_assert(std::is_same<decltype(&error_function), SUNErrHandlerFn>::value, "error_function must match SUNErrHandlerFn");


// OdeSolverCVODE definitions

OdeSolverCVODE::OdeSolverCVODE(const size_t num_dim):
    sun_context(nullptr),
    t(0.0),
    x(nullptr),
    abstol(nullptr),
    reltol(1.0e-10),
    A(nullptr),
    linear_solver(nullptr),
    cvode_mem(nullptr),
    init_error_ID(ErrorHandler::no_error),
    user_data{}
{
    // Create the SUNDIALS context
    size_t* error_ID_ptr = &(init_error_ID);
    user_data.error_ID_ptr = error_ID_ptr;    // errors are saved in this->init_error_ID
    HANDLE_ERROR_CODE(SUNContext_Create(SUN_COMM_NULL, &sun_context));
    HANDLE_ERROR_CODE(SUNContext_ClearErrHandlers(sun_context));
    HANDLE_ERROR_CODE(SUNContext_PushErrHandler(sun_context, error_function, (void*)&user_data));

    // Setup vectors
    HANDLE_RETURN_PTR(x, N_VNew_Serial(num_dim, sun_context));
    HANDLE_RETURN_PTR(abstol, N_VNew_Serial(num_dim, sun_context))
    HANDLE_RETURN_PTR(constraints, N_VNew_Serial(num_dim, sun_context));

    for (size_t i = 0; i < num_dim; i++)
        NV_Ith_S(abstol, i) = 1e-11;          // molar concentrations
    NV_Ith_S(abstol, 0) = 1e-10;              // R
    NV_Ith_S(abstol, 1) = 1e-10;              // R_dot
    NV_Ith_S(abstol, 2) = 1e-8;               // T
    NV_Ith_S(abstol, num_dim-1) = 1e-10;      // E_diss
    
    for (size_t i = 0; i < num_dim; i++)
        NV_Ith_S(constraints, i) = 1.0;      // molar concentrations c_i >= 0.0
    NV_Ith_S(constraints, 0) = 2.0;          // R > 0.0
    NV_Ith_S(constraints, 1) = 0.0;          // R_dot no constraint
    NV_Ith_S(constraints, 2) = 2.0;          // T > 0.0
    NV_Ith_S(constraints, num_dim-1) = 0.0;  // E_diss no constraint

    // Setup CVODE
    HANDLE_RETURN_PTR(cvode_mem, CVodeCreate(CV_BDF, sun_context));
    HANDLE_ERROR_CODE(CVodeInit(cvode_mem, right_hand_side, 0.0, x));
    HANDLE_ERROR_CODE(CVodeSetMaxNumSteps(cvode_mem, 10000000000));
    HANDLE_ERROR_CODE(CVodeSetMaxHnilWarns(cvode_mem, 10));    // maximum number of warnings for t+h=t
    HANDLE_ERROR_CODE(CVodeSetMaxStep(cvode_mem, 1.0e-3));     // Limit max step size to 1 ms
    HANDLE_ERROR_CODE(CVodeSetStabLimDet(cvode_mem, SUNTRUE));
    HANDLE_ERROR_CODE(CVodeSVtolerances(cvode_mem, reltol, abstol));
    HANDLE_ERROR_CODE(CVodeSetConstraints(cvode_mem, constraints));
    //HANDLE_ERROR_CODE(CVodeSStolerances(cvode_mem, 1e-10, 1e-10));

    // Setup linear solver
    HANDLE_RETURN_PTR(A, SUNDenseMatrix(num_dim, num_dim, sun_context));
    HANDLE_RETURN_PTR(linear_solver, SUNLinSol_Dense(x, A, sun_context));
    HANDLE_ERROR_CODE(CVodeSetLinearSolver(cvode_mem, linear_solver, A));
}


OdeSolverCVODE::~OdeSolverCVODE()
{
    size_t* error_ID_ptr = &(init_error_ID);
    N_VDestroy(x);
    N_VDestroy(abstol);
    N_VDestroy(constraints);
    CVodeFree(&cvode_mem);
    HANDLE_ERROR_CODE(SUNLinSolFree(linear_solver));
    SUNMatDestroy(A);
    HANDLE_ERROR_CODE(SUNContext_Free(&sun_context));
}


// used in OdeSolverCVODE::solve
void init_solve(
    void* cvode_mem,
    void* user_data,
    N_Vector x,
    size_t* error_ID_ptr
)
{
    HANDLE_ERROR_CODE(CVodeReInit(cvode_mem, 0.0, x));
    HANDLE_ERROR_CODE(CVodeSetUserData(cvode_mem, user_data));
    HANDLE_ERROR_CODE(CVodeSetInitStep(cvode_mem, 1.0e-20));
}


// used in OdeSolverCVODE::solve 
// gathers the fields (statistics) of OdeSolution from CVODE
void construct_solution(
    OdeSolution& solution,
    void* cvode_mem,
    size_t* error_ID_ptr,
    N_Vector x
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
    HANDLE_ERROR_CODE(CVodeGetNumLinIters(cvode_mem, &helper1));
    solution.num_lin_iters = helper1;
    HANDLE_ERROR_CODE(CVodeGetNumNonlinSolvIters(cvode_mem, &helper1));
    solution.num_nonlin_iters = helper1;
    HANDLE_ERROR_CODE(CVodeGetEstLocalErrors(cvode_mem, x));
    solution.total_error = std::vector<double>(NV_DATA_S(x), NV_DATA_S(x) + NV_LENGTH_S(x));
}


OdeSolution OdeSolverCVODE::solve(
    const double t_max,
    OdeFun* ode_ptr,
    double timeout,
    bool save_solution
)
{
    // Setup resources from this project
    Timer timer; timer.start();
    OdeSolution solution;
    user_data.ode_ptr = ode_ptr;
    user_data.timer_ptr = &timer;
    user_data.timeout = timeout;
    user_data.timed_out = false;
    user_data.error_ID_ptr = &(solution.error_ID);
    user_data.ID = ode_ptr->cpar.ID;

    if (init_error_ID != ErrorHandler::no_error)
    {
        solution.error_ID = init_error_ID;
        return solution;
    }
    if (ode_ptr->cpar.error_ID != ErrorHandler::no_error)
    {
        solution.error_ID = ode_ptr->cpar.error_ID;
        return solution;
    }

    // Setup CVODE resources and initial conditions
    solution.num_dim = NV_LENGTH_S(x);
    ode_ptr->initial_conditions(NV_DATA_S(x));
    solution.push_t_x(0.0, NV_DATA_S(x));
    init_solve(cvode_mem, &user_data, x, &(solution.error_ID));
    if (solution.error_ID != ErrorHandler::no_error)  return solution;

    // Solve
    int itask = save_solution ? CV_ONE_STEP : CV_NORMAL;
    while (true)
    {
        // Integration (step)
        int retval = CVode(cvode_mem, t_max, x, &t, itask);
        solution.push_t_x(t, NV_DATA_S(x));

        // Success
        if (retval == CV_SUCCESS) { }

        // Failure
        else
        {
            // Timeout and other errors handled in error_function
            if (solution.error_ID == ErrorHandler::no_error)
            {
                solution.error_ID = LOG_ERROR(Error::severity::error, Error::type::cvode,
                    "CVODE error: CVode returned with code " + std::to_string(retval) + ", but no error was logged. This is unexpected behaviour"
                );
            }
        }

        // Exit conditions
        if (t >= t_max)  break;
        if (user_data.timed_out)  break;
        if (retval != CV_SUCCESS)  break;
    }

    // fill solution
    construct_solution(solution, cvode_mem, &(solution.error_ID), x);
    solution.runtime       = timer.lap();
    user_data.ode_ptr      = nullptr;
    user_data.timer_ptr    = nullptr;
    user_data.error_ID_ptr = &init_error_ID;

    return solution;
}
#include "ode_solver_sundials.h"

#include <type_traits>
static_assert(std::is_same<sunrealtype, double>::value, "sunrealtype must be double");


#define HANDLE_ERROR_CODE(...) \
{ \
    int retval = __VA_ARGS__; \
    if (retval != 0) \
    { \
        this->error_ID = LOG_ERROR(Error::severity::error, Error::type::preprocess, "CVODE error: " # __VA_ARGS__ "returned with code " + std::to_string(retval)); \
        return; \
    } \
}

#define HANDLE_RETURN_PTR(pointer, ...) \
{ \
    pointer = __VA_ARGS__; \
    if (pointer == nullptr) \
    { \
        this->error_ID = LOG_ERROR(Error::severity::error, Error::type::preprocess, "CVODE error: " # __VA_ARGS__ "returned with nullptr"); \
        return; \
    } \
}


int right_hand_side(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
    OdeFun* ode = (OdeFun*)user_data;
    is_success ret = ode->operator()(t, NV_DATA_S(y), NV_DATA_S(ydot));
    return !ret;
}


OdeSolverCVODE::OdeSolverCVODE(OdeFun *ode):
    sun_context(nullptr),
    t(0.0),
    x(nullptr),
    abstol(nullptr),
    A(nullptr),
    linear_solver(nullptr),
    cvode_mem(nullptr)
{
    const double _abstol = 1e-10;
    const double _reltol = 1e-10;
    size_t num_dim = ode->par->num_species+4;

    // Create the SUNDIALS context
    HANDLE_ERROR_CODE(SUNContext_Create(SUN_COMM_NULL, &sun_context));

    // Setup vectors
    HANDLE_RETURN_PTR(x, N_VNew_Serial(num_dim, sun_context));
    ode->initial_conditions(NV_DATA_S(x));
    HANDLE_RETURN_PTR(abstol, N_VNew_Serial(num_dim, sun_context))
    for (int i = 0; i < num_dim; i++)
        NV_Ith_S(abstol, i) = _abstol;

    // Setup CVODE
    HANDLE_RETURN_PTR(cvode_mem, CVodeCreate(CV_BDF, sun_context));
    HANDLE_ERROR_CODE(CVodeInit(cvode_mem, right_hand_side, 0.0, x));
    HANDLE_ERROR_CODE(CVodeSetUserData(cvode_mem, ode));
    HANDLE_ERROR_CODE(CVodeSetMaxNumSteps(cvode_mem, 10000000));
    HANDLE_ERROR_CODE(CVodeSVtolerances(cvode_mem, _reltol, abstol));

    // Setup linear solver
    HANDLE_RETURN_PTR(A, SUNDenseMatrix(num_dim, num_dim, sun_context));
    HANDLE_RETURN_PTR(linear_solver, SUNLinSol_Dense(x, A, sun_context));
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
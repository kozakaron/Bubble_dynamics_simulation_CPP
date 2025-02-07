#ifdef BENCHMARK

#include "common.h"
#include "parameters.h"
#include "control_parameters.h"
#include "ode_fun.h"
#include "ode_solver.h"
#include "test_list.h"

#include <thread>
#include <vector>
#include <iostream>
#include <iomanip>
#include <atomic>

constexpr size_t num_tasks = 64;
std::atomic<size_t> task_counter;

void task()
{
    // Runtime of solver ~ 2.3 sec
    // Runtime of setup ~ 0.000002 sec (2 us)
    ControlParameters cpar;
    OdeFun ode;
    double x[200];
    ode.init(cpar);
    ode.initial_conditions(x);
    auto ode_fun = [&ode](const double t, const double *x, double *dxdt) -> is_success { return ode(t, x, dxdt); };
    RKCK45 solver;

    while(true)
    {
        const size_t task_ID = task_counter.fetch_add(1, std::memory_order_relaxed);
        if (task_ID >= num_tasks) break;

        solver.solve(0.0, 100.0e-6, (double*)x, ode.par->num_species+4, ode_fun, &ode.error_ID, 60.0, false);
        //OdeSolution sol = solver.get_solution();
        //std::cout << task_ID << ". " << sol.runtime << "\n";
    }
}

void benchmark_speedup()
{
    std::cout << colors::bold << "Measure multithreaded runtime with RKCK45 solver and simple chemkin_ar_he mechanism" << colors::reset << std::endl;
    std::cout << "    Number of threads: " << std::thread::hardware_concurrency() << std::endl;
    std::cout << "    Single-threaded runtime of one task: ";

    Timer timer;
    timer.start();
    task_counter = num_tasks - 1;
    task();
    double single_thread_runtime = timer.lap();

    std::cout << single_thread_runtime << " s;   Number of tasks: " << num_tasks << std::endl;
    std::cout << "    Threads [-]   |   Runtime [s]   |   ~Speedup [-]\n";
    single_thread_runtime *= num_tasks;

    for (size_t num_threads = std::thread::hardware_concurrency(); num_threads >= 1; num_threads /= 2)
    {
        {
            std::vector<std::jthread> threads(num_threads);
            task_counter = 0;

            timer.start();

            for (size_t i = 0; i < num_threads; i++)
            {
                threads[i] = std::jthread(task);
            }
        }

        double runtime = timer.lap();
        std::cout << "    " << std::setw(8) << num_threads << "      |   " << std::setw(8) << runtime;
        std::cout  << "      |   " << single_thread_runtime / runtime << std::endl;
    }

}

#endif // BENCHMARK
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

size_t num_tasks;
std::atomic<size_t> task_counter;
constexpr size_t task_per_thread = 4;

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
    std::cout << "    Threads [-]     |  Runtime/task [s] |     Speedup [-]\n";

    Timer timer;
    double single_thread_runtime;
    for (size_t num_threads = 1; num_threads <= std::thread::hardware_concurrency(); num_threads *= 2)
    {
        {
            std::vector<std::jthread> threads(num_threads);
            task_counter = 0;
            num_tasks = task_per_thread * num_threads;

            timer.start();

            for (size_t i = 0; i < num_threads; i++)
            {
                threads[i] = std::jthread(task);
            }
        }

        double runtime = timer.lap() / num_threads / task_per_thread;
        if (num_threads == 1) single_thread_runtime = runtime;
        std::cout << "    " << std::setw(10) << num_threads << "      |   " << std::setw(10) << runtime;
        std::cout  << "      |     " << single_thread_runtime / runtime << std::endl;
    }
}

#endif // BENCHMARK
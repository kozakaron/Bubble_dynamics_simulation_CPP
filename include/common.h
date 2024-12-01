#ifndef COMMON_H
#define COMMON_H
#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <array>
#include <cmath>

#define WHITE "\033[0m"
#define RED "\033[31m"
#define GREEN "\033[32m"
#define YELLOW "\033[33m"
#define BLUE "\033[94m"
#define MAGENTA "\033[35m"
#define CYAN "\033[36m"
#define BOLD "\033[1m"
#define ITALIC "\033[3m"
#define UNDERLINE "\033[4m"
#define REVERSED "\033[7m"
#define RESET "\033[0m"

// cout overload for 1D std::array
template <typename T, size_t N>
std::ostream &operator<<(std::ostream &os, const std::array<T, N> &arr);

// cout overload for 2D std::array
template <typename T, size_t N, size_t M>
std::ostream &operator<<(std::ostream &os, const std::array<std::array<T, M>, N> &arr);

class Timer
{
public:
    Timer() {}
    // Start the timer
    void start();
    // Measure the time elapsed since the last start
    double lap();

private:
    std::chrono::high_resolution_clock::time_point start_time;
};

#ifdef BENCHMARK
#define BENCHMARK_LINE(line, N) \
    { \
        Timer benchmark_timer; \
        size_t run = 0; \
        benchmark_timer.start(); \
        for (; run < N; ++run) \
        { \
            line; \
        } \
        double benchmark_time = benchmark_timer.lap(); \
        std::cout << "Runtime of " << BLUE << #line << RESET; \
        std::cout << ":  " << BOLD << 1e6 * benchmark_time / N << " us" << RESET << std::endl; \
    }
#else
#define BENCHMARK_LINE(line, N)
#endif

#endif // COMMON_H
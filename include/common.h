#ifndef COMMON_H
#define COMMON_H
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <chrono>
#include <array>
#include <vector>
#include <algorithm>
#include <cmath>
#include <mutex>

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

// Timer class for benchmarking
class Timer
{
public:
    Timer() {}
    // Start the timer
    void start();
    // Measure the time elapsed since the last start in seconds
    double lap();
    // Get current time
    static std::string current_time();

private:
    std::chrono::high_resolution_clock::time_point start_time;
};

class Error
{
public:
    std::string message;
    std::string function;
    std::string file;
    std::string time;
    size_t line;
    size_t ID;

    Error();
    Error(
        const std::string &message,
        const std::string &function,
        const std::string &file,
        const size_t line,
        const size_t ID = 0
    );
    std::string to_string(const bool color = false) const;
    friend std::ostream &operator<<(std::ostream &os, const Error &err);
};

class ErrorHandler
{
private:
    static std::vector<Error> errors;
    static std::mutex mutex;
    static std::ofstream log_file;
public:
    static bool print_when_log;

    static void set_log_file(const std::string &filename);
    static void log_error(const Error &err);
    static void print_errors();
    static size_t get_error_count();
    static void clear_errors();
    static Error get_error(size_t index);
};

// Error macro: ERROR(string message, size_t ID=0)
#define ERROR(message, ID) Error(message, __FUNCTION__, __FILE__, __LINE__, ID)
#define LOG_ERROR(message, ID) ErrorHandler::log_error(ERROR(message, ID))

// Benchmarking macro
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
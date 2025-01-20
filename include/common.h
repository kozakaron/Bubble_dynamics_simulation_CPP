#ifndef COMMON_H
#define COMMON_H
#include <iostream>
#include <string>
#include <chrono>
#include <vector>
#include <mutex>

typedef bool is_success;

namespace colors
{
    const std::string white =       "\033[0m";
    const std::string red =         "\033[31m";
    const std::string green =       "\033[32m";
    const std::string yellow =      "\033[33m";
    const std::string blue =        "\033[94m";
    const std::string magenta =     "\033[35m";
    const std::string cyan =        "\033[36m";
    const std::string bold =        "\033[1m";
    const std::string italic =      "\033[3m";
    const std::string underline =   "\033[4m";
    const std::string reversed =    "\033[7m";
    const std::string reset =       "\033[0m";
}

template <typename T>
std::string to_string(T* arr, size_t len);

template <typename T>
std::string to_string(T* arr, size_t len1, size_t len2);

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
    // Format elapsed time (0.002 -> "2 ms", 0.000000002 -> "2 ns")
    static std::string format_time(double time);

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
    static constexpr size_t no_error = std::numeric_limits<size_t>::max();

    static void set_log_file(const std::string &filename);
    static size_t log_error(const Error &err);
    static void print_errors();
    static size_t get_error_count();
    static void clear_errors();
    static Error get_error(size_t index);
};

// Error macro: ERROR(string message, size_t ID=0)
#define ERROR(message, ID) Error(message, __FUNCTION__, __FILE__, __LINE__, ID)
#define LOG_ERROR(message, ID) ErrorHandler::log_error(ERROR(message, ID))

// Unrolling macro, optimization off macro
#if defined(__GNUC__) && !defined(__clang__)
#define LOOP_UNROLL(n) _Pragma("GCC unroll " #n)
#define NO_OPTIMIZATION _Pragma("GCC optimize(\"O0\")")
#elif defined(__clang__)
#define LOOP_UNROLL(n) _Pragma("clang loop unroll_count(" #n ")")
#define NO_OPTIMIZATION _Pragma("clang optimize off")
#else
#define LOOP_UNROLL(n)
#define NO_OPTIMIZATION
#endif

// Benchmarking macro
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
        double run_time = benchmark_time / N; \
        std::cout << "Runtime of " << colors::blue << #line << colors::reset << " is " << colors::bold << Timer::format_time(run_time) << colors::reset << std::endl; \
    }

#endif // COMMON_H
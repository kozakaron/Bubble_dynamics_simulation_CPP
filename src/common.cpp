#include <fstream>
#include <algorithm>
#include <sstream>
#include <type_traits>
#include <iomanip>

#include "common.h"

template <typename T>
std::string to_string(T* arr, size_t len)
{
    std::stringstream ss;
    ss << "{";
    using ElementType = typename std::decay<decltype(arr[0])>::type;
    bool is_float = std::is_floating_point<ElementType>::value;
    for (size_t i = 0; i < len; ++i)
    {
        if (is_float)
            ss << std::scientific << std::setprecision(std::numeric_limits<ElementType>::max_digits10);
        ss << arr[i];
        if (i < len - 1)
            ss << ", ";
    }
    ss << "}";
    return ss.str();
}

template std::string to_string(float* arr, size_t len);
template std::string to_string(double* arr, size_t len);
template std::string to_string(char* arr, size_t len);
template std::string to_string(short* arr, size_t len);
template std::string to_string(int* arr, size_t len);
template std::string to_string(long* arr, size_t len);
template std::string to_string(long long* arr, size_t len);
template std::string to_string(unsigned short* arr, size_t len);
template std::string to_string(unsigned* arr, size_t len);
template std::string to_string(unsigned long* arr, size_t len);
template std::string to_string(unsigned long long* arr, size_t len);
template std::string to_string(std::string* arr, size_t len);


template <typename T>
std::string to_string(T** arr, size_t len1, size_t len2)
{
    std::stringstream ss;
    ss << "{\n";
    for(size_t i = 0; i < len1; ++i)
    {
        ss << "    " << to_string(arr[i], len2);
        if (i < len1 - 1)
            ss << ",\n";
        else
            ss << "\n";
    }
    ss << "}";
    return ss.str();
}

template std::string to_string(float** arr, size_t len1, size_t len2);
template std::string to_string(double** arr, size_t len1, size_t len2);
template std::string to_string(char** arr, size_t len1, size_t len2);
template std::string to_string(short** arr, size_t len1, size_t len2);
template std::string to_string(int** arr, size_t len1, size_t len2);
template std::string to_string(long** arr, size_t len1, size_t len2);
template std::string to_string(long long** arr, size_t len1, size_t len2);
template std::string to_string(unsigned short** arr, size_t len1, size_t len2);
template std::string to_string(unsigned** arr, size_t len1, size_t len2);
template std::string to_string(unsigned long** arr, size_t len1, size_t len2);
template std::string to_string(unsigned long long** arr, size_t len1, size_t len2);
template std::string to_string(std::string** arr, size_t len1, size_t len2);


void Timer::start()
{
    start_time = std::chrono::high_resolution_clock::now();
}

double Timer::lap()
{
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
    return 1e-9 * static_cast<double>(duration);
}

std::string Timer::current_time()
{
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    std::tm *buf;
    std::stringstream ss;
#ifdef _WIN32
    buf = new std::tm;
    localtime_s(buf, &in_time_t);
#else
    buf = localtime(&in_time_t);
#endif
    ss << std::put_time(buf, "%Y.%m.%d %X");
    return ss.str();
}

std::string Timer::format_time(double time)
{
    std::stringstream ss;
    if (time < 1e-6)
    {
        ss << std::setprecision(2) << std::fixed << time * 1e9 << " ns";
    }
    else if (time < 1e-3)
    {
        ss << std::setprecision(2) << std::fixed << time * 1e6 << " us";
    }
    else if (time < 1)
    {
        ss << std::setprecision(2) << std::fixed << time * 1e3 << " ms";
    }
    else
    {
        ss << std::setprecision(3) << std::fixed << time << " s";
    }
    return ss.str();
}

Error::Error()
{
    this->message = "";
    this->function = "";
    this->file = "";
    this->line = 0;
    this->ID = 0;
    this->time = "";
}

Error::Error(
    const std::string &message,
    const std::string &function,
    const std::string &file,
    const size_t line,
    const size_t ID
)
{
    this->message = message;
    this->function = function;
    this->file = file;
    std::replace(this->file.begin(), this->file.end(), '\\', '/');
    this->line = line;
    this->ID = ID;
    this->time = Timer::current_time();
}

std::string Error::to_string(const bool color) const
{
    std::stringstream ss;
    ss << this->time << ": ";
    if (color) ss << colors::red;
    ss << this->message;
    if (color) ss << colors::reset;
    ss << " in " << this->file << ":" << this->line << ": " << this->function << "();";
    if (this->ID != 0) ss << " ID=" << this->ID << ";";
    return ss.str();
}

std::ostream &operator<<(std::ostream &os, const Error &err)
{
    os << err.to_string(true);
    return os;
}

std::vector<Error> ErrorHandler::errors;
std::mutex ErrorHandler::mutex;
std::ofstream ErrorHandler::log_file;
bool ErrorHandler::print_when_log = true;

void ErrorHandler::set_log_file(const std::string &filename)
{
    ErrorHandler::log_file.open(filename, std::ios::out | std::ios::app);
    if (!ErrorHandler::log_file.is_open())
    {
        LOG_ERROR("Could not open log file " + filename, 0);
    }
}

size_t ErrorHandler::log_error(const Error &err)
{
    if (ErrorHandler::print_when_log)
    {
        std::cerr << err << std::endl;
    }
    if (ErrorHandler::log_file.is_open())
    {
        ErrorHandler::log_file << err.to_string(false) << std::endl;
    }
    std::lock_guard<std::mutex> lock(ErrorHandler::mutex);
    ErrorHandler::errors.push_back(err);
    return ErrorHandler::errors.size() - 1;
}

void ErrorHandler::print_errors()
{
    std::lock_guard<std::mutex> lock(ErrorHandler::mutex);
    if (ErrorHandler::errors.empty())
    {
        std::cout << colors::bold << colors::green << "No errors." << colors::reset << std::endl;
        return;
    }
    std::cout << colors::bold << colors::red << "Errors:" << colors::reset << std::endl;
    for (const auto &err : errors)
    {
        std::cout << "    " << err.to_string(true) << std::endl;
    }
}

size_t ErrorHandler::get_error_count()
{
    std::lock_guard<std::mutex> lock(ErrorHandler::mutex);
    return ErrorHandler::errors.size();
}

void ErrorHandler::clear_errors()
{
#ifndef TEST
    std::cerr << colors::bold << colors::red << "It is not recommended to clear errors outside test mode. " << colors::reset << std::endl;
    // Why? becouse error_ID is used to point to certain errors
#endif
    std::lock_guard<std::mutex> lock(ErrorHandler::mutex);
    ErrorHandler::errors.clear();
}

Error ErrorHandler::get_error(size_t index)
{
    if (index == ErrorHandler::no_error)
        return ERROR("No error occured", 0);
    if (index >= ErrorHandler::errors.size())
        return ERROR("Error index " + std::to_string(index) + " out of bounds, number of errors: " + std::to_string(ErrorHandler::errors.size()) +\
                     " (Perhaps you used ErrorHandler::clear_errors)", 0);
    std::lock_guard<std::mutex> lock(ErrorHandler::mutex);
    return ErrorHandler::errors[index];
}
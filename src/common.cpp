#include "common.h"

template <typename T, size_t N>
std::ostream  &operator<<(std::ostream  &os, const std::array<T, N> &arr)
{
    os << std::setfill(' ') << "[";
    for (size_t i = 0; i < N; ++i)
    {
        os << std::setw(8) << arr[i];
        if (i < N - 1)
        {
            os << ",";
        }
    }
    os << "]";
    return os;
}


template <typename T, size_t N, size_t M>
std::ostream  &operator<<(std::ostream  &os, const std::array<std::array<T, M>, N> &arr)
{
    os << "[";
    for (size_t i = 0; i < N; ++i)
    {
        os << "\n  ";
        os << arr[i];
        if (i < N - 1)
        {
            os << ", ";
        }
    }
    os << "\n]";
    return os;
}

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
    if (color) ss << RED;
    ss << this->message;
    if (color) ss << RESET;
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
        LOG_ERROR("Could not open log file " + filename, 1);
    }
}

void ErrorHandler::log_error(const Error &err)
{
    std::lock_guard<std::mutex> lock(ErrorHandler::mutex);
    ErrorHandler::errors.push_back(err);
    if (ErrorHandler::print_when_log)
    {
        std::cerr << err << std::endl;
    }
    if (ErrorHandler::log_file.is_open())
    {
        ErrorHandler::log_file << err.to_string(false) << std::endl;
    }
}

void ErrorHandler::print_errors()
{
    std::lock_guard<std::mutex> lock(ErrorHandler::mutex);
    if (ErrorHandler::errors.empty())
    {
        std::cout << BOLD << GREEN << "No errors." << RESET << std::endl;
        return;
    }
    std::cout << BOLD << RED << "Errors:" << RESET << std::endl;
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
    std::lock_guard<std::mutex> lock(ErrorHandler::mutex);
    ErrorHandler::errors.clear();
}

Error ErrorHandler::get_error(size_t index)
{
    std::lock_guard<std::mutex> lock(ErrorHandler::mutex);
    return ErrorHandler::errors[index];
}
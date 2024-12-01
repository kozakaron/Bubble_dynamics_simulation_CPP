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
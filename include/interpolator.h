#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include <vector>
#include <string>
#include <tuple>

#include "common.h"


// Class used to load CSV data and interpolate continous
//    * R(t) profile instead of built-in bubble dynamics
//    * p(t) profile instead of built-in excitations
// Can also provide smooth first and second derivatives
class Interpolator {
public:
    std::vector<double> t_data;
    std::vector<double> x_data;
	std::vector<double> v_data; // R_dot
    size_t error_ID;
    std::string filename;  // Stored filename with normalized path separators (/ only)

    // Construct from CSV path with the following features:
    //    1. Remove not strictly monotanious points
    //    2. Remove start-up points with exponentially growing steps
    // Expected format for 1:
    //     t [s], R [m]
    //     0.0,1e-5
    //     ...
    // Expected format for 2:
    //     t [s], p [Pa]
    //     0.0,1e5
    //     ...
    Interpolator(std::string csv_path);
    Interpolator();

    // Interpolates between read data points. Returns interpolated value and first two derivatives.
    // Uses quartic (4th order) Largrange polynomial interpolation using 5 data points.
    std::tuple<double, double, double> interpolate(const double t);
};



#endif // INTERPOLATOR_H
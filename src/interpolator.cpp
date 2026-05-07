#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <span>
#include <filesystem>

#include "common.h"
#include "interpolator.h"

/*
// Helper function to print expected CSV format
static void print_interpolator_format_help()
{
    std::cout << colors::bold << "Expected CSV format:" << colors::reset << "\n"
              << " * Header: 1 line (e.g., 't [s], R [m]' or 't [s], p [Pa]')\n"
              << " * Seperator: comma (',')\n"
              << " * First column: time values in seconds (strictly monotonically increasing)\n"
              << " * Second column: corresponding data values (bubble radius in meters or ambient pressure in pascals)\n"
              << " * At least 10 data points required\n"
              << " * Example:\n"
              << "      t [s],R [m]\n"
              << "      0.0,1.0e-5\n"
              << "      1.0e-6,1.001e-5\n"
              << "      2.0e-6,1.002e-5\n"
              << "      ...\n" << std::endl;
}


Interpolator::Interpolator()
{
    this->error_ID = ErrorHandler::no_error;
}


Interpolator::Interpolator(std::string csv_path)
{
    Timer timer; timer.start();
    this->error_ID = ErrorHandler::no_error;
    
    // Convert to absolute path
    std::filesystem::path abs_path = std::filesystem::absolute(csv_path);
    this->filename = abs_path.string();
    
    // Replace \ and \\ with / for consistent path representation
    for (size_t i = 0; i < this->filename.length(); ++i)
    {
        if (this->filename[i] == '\\')
        {
            this->filename[i] = '/';
        }
    }
    size_t pos = 0;
    while ((pos = this->filename.find("//", pos)) != std::string::npos)
    {
        this->filename.replace(pos, 2, "/");
    }
    
    // Check if file exists
    std::ifstream file(csv_path);
    if (!file.is_open())
    {
        this->error_ID = LOG_ERROR("Failed to open CSV file: " + csv_path);
        return;
    }
    
    std::string line;
    std::vector<double> t_temp, x_temp;
    
    // Skip header row
    if (!std::getline(file, line))
    {
        print_interpolator_format_help();
        this->error_ID = LOG_ERROR("CSV file is empty: " + this->filename);
        file.close();
        return;
    }
    
    // Read data rows
    size_t row_count = 2;
    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string col1, col2;
        
        // Read columns
        if (!std::getline(ss, col1, ',') || !std::getline(ss, col2))
        {
            print_interpolator_format_help();
            this->error_ID = LOG_ERROR(
                "Not enough columns in row " + std::to_string(row_count) + " (\"" + line + "\") in \"" + this->filename + "\""
            );
            file.close();
            return;
        }
        
        // Parse doubles
        double t_val, x_val;
        try
        {
            t_val = std::stod(col1);
            x_val = std::stod(col2);
        }
        catch (const std::exception& e)
        {
            print_interpolator_format_help();
            this->error_ID = LOG_ERROR(
                "Failed to parse doubles in row " + std::to_string(row_count) + " (\"" + line + "\") in \"" + this->filename + "\""
            );
            file.close();
            return;
        }
        
        t_temp.push_back(t_val);
        x_temp.push_back(x_val);
        row_count++;
    }
    file.close();
    
    // Check minimum data points
    if (t_temp.size() < 10)
    {
        print_interpolator_format_help();
        this->error_ID = LOG_ERROR(
            "At least 10 data points are required, but only got " + std::to_string(t_temp.size()) + " in \"" + this->filename + "\""
        );
        return;
    }

    // Remove exponentially growing initial steps, if present
    const double max_spacing_ratio = 50.0;
    double x_min = x_temp[0], x_max = x_temp[0];
    for (size_t i = 1; i < x_temp.size(); ++i)
    {
        if (x_temp[i] < x_min) x_min = x_temp[i];
        if (x_temp[i] > x_max) x_max = x_temp[i];
    }
    const double data_range = std::abs(x_max - x_min);
    const double quiet_tol = 1e-8 * data_range;
    
    size_t first_good = 0;
    for (size_t i = 0; i < t_temp.size() - 4; ++i)
    {
        // Calculate 4 gaps for a 5-point stencil
        double dt_min = t_temp[i + 1] - t_temp[i];
        double dt_max = dt_min;
        for (size_t j = i + 1; j < i + 4; ++j)
        {
            double dt = t_temp[j + 1] - t_temp[j];
            if (dt < dt_min) dt_min = dt;
            if (dt > dt_max) dt_max = dt;
        }
        
        const double ratio = (dt_min > 0) ? (dt_max / dt_min) : 1e10;
        const double val_change = std::abs(x_temp[i + 4] - x_temp[i]);
        
        if (ratio > max_spacing_ratio && val_change < quiet_tol)
        {
            first_good = i + 1;  // drop this point
        }
        else
        {
            break;
        }
    }
    
    if (first_good > 0)
    {
        LOG_ERROR(
            Error::severity::info,
            Error::type::preprocess,
            "Removed " + std::to_string(first_good) + " ill-conditioned startup points from \"" + this->filename + "\"" 
        );
        
        // Shift data to remove startup points
        std::vector<double> t_filtered(t_temp.begin() + first_good, t_temp.end());
        std::vector<double> x_filtered(x_temp.begin() + first_good, x_temp.end());
        t_temp = std::move(t_filtered);
        x_temp = std::move(x_filtered);
    }
    
    // Remove non-monotonous points (strict inequality to prevent duplicates)
    size_t removed_count = 0;
    this->t_data.reserve(t_temp.size());
    this->x_data.reserve(t_temp.size());
    this->t_data.push_back(t_temp[0]);
    this->x_data.push_back(x_temp[0]);
    
    for (size_t i = 1; i < t_temp.size(); ++i)
    {
        if (t_temp[i] > this->t_data.back())  // Changed from >= to >
        {
            this->t_data.push_back(t_temp[i]);
            this->x_data.push_back(x_temp[i]);
        }
        else
        {
            removed_count++;
        }
    }
    
    if (removed_count > 0)
    {
        LOG_ERROR(
            Error::severity::info,
            Error::type::preprocess,
            "Removed " + std::to_string(removed_count) + " not strictly monotonious points from \"" + this->filename + "\"" 
        );
    }
    
    // Log success with timing
    const double runtime = timer.lap();
    LOG_ERROR(
        Error::severity::info,
        Error::type::preprocess,
        "Loaded " + std::to_string(this->t_data.size()) + " data points from \"" + this->filename + "\" in " + Timer::format_time(runtime)
    );
}


// Lagrange basis polynomial L_i(t) for i-th point among 5 stencil points
inline double lagrange_basis_L(const double t, const double ti, const double tj, const double tk, const double tl, const double tm)
{
    return (t - tj) * (t - tk) * (t - tl) * (t - tm) /
           ((ti - tj) * (ti - tk) * (ti - tl) * (ti - tm));
}


// First derivative of Lagrange basis polynomial L_i'(t)
inline double lagrange_basis_dL(const double t, const double ti, const double tj, const double tk, const double tl, const double tm)
{
    const double denom = (ti - tj) * (ti - tk) * (ti - tl) * (ti - tm);
    return ((t - tk) * (t - tl) * (t - tm) +
            (t - tj) * (t - tl) * (t - tm) +
            (t - tj) * (t - tk) * (t - tm) +
            (t - tj) * (t - tk) * (t - tl)) / denom;
}


// Second derivative of Lagrange basis polynomial L_i''(t)
inline double lagrange_basis_ddL(const double t, const double ti, const double tj, const double tk, const double tl, const double tm)
{
    const double denom = (ti - tj) * (ti - tk) * (ti - tl) * (ti - tm);
    return 2.0 * ((t - tl) * (t - tm) + (t - tk) * (t - tm) + (t - tk) * (t - tl) +
                  (t - tj) * (t - tm) + (t - tj) * (t - tl) + (t - tj) * (t - tk)) / denom;
}


// Quartic Lagrange polynomial evaluation at point t using 5-point stencil
// Returns interpolated value P(t), first derivative P'(t), and second derivative P''(t)
// Uses explicit formulation of Lagrange basis polynomials for clarity and performance
// Implements local time normalization to [0,1] for improved numerical conditioning
inline std::tuple<double, double, double> lagrange_quartic_raw(
    const double t,
    const std::span<const double> t_pts,
    const std::span<const double> x_pts
)
{
    const double t0 = t_pts[0], t1 = t_pts[1], t2 = t_pts[2], t3 = t_pts[3], t4 = t_pts[4];
    const double y0 = x_pts[0], y1 = x_pts[1], y2 = x_pts[2], y3 = x_pts[3], y4 = x_pts[4];
    
    // Local time normalization: map [t_min, t_max] → [0, 1]
    // This dramatically improves numerical conditioning when timesteps vary by orders of magnitude
    const double t_min = t0;  // First point in stencil
    const double t_max = t4;  // Last point in stencil
    const double dt_range = t_max - t_min;
    
    // Normalize query point and stencil points to [0, 1]
    const double t_norm = (t - t_min) / dt_range;
    const double t0_norm = 0.0;  // (t0 - t_min) / dt_range = 0
    const double t4_norm = 1.0;  // (t4 - t_min) / dt_range = 1
    const double t1_norm = (t1 - t_min) / dt_range;
    const double t2_norm = (t2 - t_min) / dt_range;
    const double t3_norm = (t3 - t_min) / dt_range;
    
    // Interpolated value: P(t) = Σ y_i * L_i(t)
    // Value is scale-invariant, so no adjustment needed
    const double d0  = y0 * lagrange_basis_L(t_norm, t0_norm, t1_norm, t2_norm, t3_norm, t4_norm) + 
                       y1 * lagrange_basis_L(t_norm, t1_norm, t0_norm, t2_norm, t3_norm, t4_norm) +
                       y2 * lagrange_basis_L(t_norm, t2_norm, t0_norm, t1_norm, t3_norm, t4_norm) + 
                       y3 * lagrange_basis_L(t_norm, t3_norm, t0_norm, t1_norm, t2_norm, t4_norm) +
                       y4 * lagrange_basis_L(t_norm, t4_norm, t0_norm, t1_norm, t2_norm, t3_norm);
    
    // First derivative: P'(t) = Σ y_i * L_i'(t)
    // Chain rule: dP/dt_actual = dP/dt_norm × dt_norm/dt_actual = dP/dt_norm / dt_range
    const double d1_norm = y0 * lagrange_basis_dL(t_norm, t0_norm, t1_norm, t2_norm, t3_norm, t4_norm) + 
                           y1 * lagrange_basis_dL(t_norm, t1_norm, t0_norm, t2_norm, t3_norm, t4_norm) +
                           y2 * lagrange_basis_dL(t_norm, t2_norm, t0_norm, t1_norm, t3_norm, t4_norm) + 
                           y3 * lagrange_basis_dL(t_norm, t3_norm, t0_norm, t1_norm, t2_norm, t4_norm) +
                           y4 * lagrange_basis_dL(t_norm, t4_norm, t0_norm, t1_norm, t2_norm, t3_norm);
    const double d1 = d1_norm / dt_range;
    
    // Second derivative: P''(t) = Σ y_i * L_i''(t)
    // Chain rule: d²P/dt_actual² = d²P/dt_norm² × (dt_norm/dt_actual)² = d²P/dt_norm² / dt_range²
    const double d2_norm = y0 * lagrange_basis_ddL(t_norm, t0_norm, t1_norm, t2_norm, t3_norm, t4_norm) + 
                           y1 * lagrange_basis_ddL(t_norm, t1_norm, t0_norm, t2_norm, t3_norm, t4_norm) +
                           y2 * lagrange_basis_ddL(t_norm, t2_norm, t0_norm, t1_norm, t3_norm, t4_norm) + 
                           y3 * lagrange_basis_ddL(t_norm, t3_norm, t0_norm, t1_norm, t2_norm, t4_norm) +
                           y4 * lagrange_basis_ddL(t_norm, t4_norm, t0_norm, t1_norm, t2_norm, t3_norm);
    const double d2 = d2_norm / (dt_range * dt_range);
    
    return std::make_tuple(d0, d1, d2);
}


// Blended quartic Lagrange interpolation with smooth C^2 stencil transitions
// Addresses the stencil discontinuity problem: when the interpolation point moves from one
// 5-point stencil to another, abrupt changes can occur. This function smoothly blends two
// overlapping stencils using Perlin smootherstep weighting (5th-order polynomial) to ensure
// continuity in value, first derivative, and second derivative across stencil boundaries.
// The transition occurs in the overlap region between adjacent stencils.
std::tuple<double, double, double> lagrange_quartic_blended(
    const double t,
    const std::vector<double>& t_data,
    const std::vector<double>& x_data
)
{
    const size_t n = t_data.size();
    
    // Binary search to locate interval: t[j] <= t < t[j+1]
    size_t low = 0, high = n - 1;
    while (low <= high)
    {
        size_t mid = low + (high - low) / 2;
        if (t_data[mid] <= t)
        {
            low = mid + 1;
        }
        else
        {
            high = mid - 1;
        }
    }
    size_t j = high;
    j = std::min(std::max(j, static_cast<size_t>(0)), n - 2);
    
    // Primary stencil i0, centered around j
    size_t i0 = std::min(std::max((int)j - 2, 0), (int)(n - 5));
    
    // Shifted stencil i0b (one step to the right)
    size_t i0b = std::min(i0 + 1, n - 5);
    
    // EDGE CASE 1: Stencils are identical (happens near right boundary)
    // When i0 >= n-4, there's no room for a second stencil: i0b = min(i0+1, n-5) = n-5 = i0
    // Example: n=10 → max stencil start is n-5=5. If i0=5, then i0b=min(6,5)=5 → same stencil
    // This happens when we're in the last 5 points of data. No blending possible, use single stencil.
    if (i0 == i0b)
    {
        std::span<const double> t_pts(t_data.data() + i0, 5);
        std::span<const double> x_pts(x_data.data() + i0, 5);
        return lagrange_quartic_raw(t, t_pts, x_pts);
    }
    
    // Define blending transition interval [t_lo, t_hi]
    // Stencil A: [i0, i0+1, i0+2, i0+3, i0+4]
    // Stencil B: [i0b, i0b+1, i0b+2, i0b+3, i0b+4] = [i0+1, i0+2, i0+3, i0+4, i0+5]
    // Overlap points: i0+2, i0+3, i0+4 (shared by both stencils)
    const double t_lo = t_data[i0 + 2];  // Start of transition zone
    const double t_hi = t_data[i0 + 3];  // End of transition zone
    
    // Normalize position within transition zone to s ∈ [0, 1]
    double s = (t - t_lo) / (t_hi - t_lo);
    s = std::min(std::max(s, 0.0), 1.0);
    
    // EDGE CASE 2: Query point at or before transition start (s = 0)
    // When t ≤ t_lo, we're entirely in stencil A's domain before the blend zone
    // Blending weight w would be 0 anyway (since w(0)=0), so skip computation and use stencil A
    if (s <= 0.0)
    {
        std::span<const double> t_pts(t_data.data() + i0, 5);
        std::span<const double> x_pts(x_data.data() + i0, 5);
        return lagrange_quartic_raw(t, t_pts, x_pts);
    }
    
    // EDGE CASE 3: Query point at or after transition end (s = 1)
    // When t ≥ t_hi, we're entirely in stencil B's domain after the blend zone
    // Blending weight w would be 1 anyway (since w(1)=1), so skip computation and use stencil B
    if (s >= 1.0)
    {
        std::span<const double> t_pts(t_data.data() + i0b, 5);
        std::span<const double> x_pts(x_data.data() + i0b, 5);
        return lagrange_quartic_raw(t, t_pts, x_pts);
    }
    
    // Perlin smootherstep: C^2 smooth weight function
    // w = s^3(6s^2 - 15s + 10)
    const double w = s * s * s * (s * (s * 6.0 - 15.0) + 10.0);
    
    // Evaluate both stencils
    std::span<const double> t_pts0(t_data.data() + i0, 5);
    std::span<const double> x_pts0(x_data.data() + i0, 5);
    auto [d0_0, d1_0, d2_0] = lagrange_quartic_raw(t, t_pts0, x_pts0);
    
    std::span<const double> t_pts1(t_data.data() + i0b, 5);
    std::span<const double> x_pts1(x_data.data() + i0b, 5);
    auto [d0_1, d1_1, d2_1] = lagrange_quartic_raw(t, t_pts1, x_pts1);
    
    // Blend results
    const double d0 = (1.0 - w) * d0_0 + w * d0_1;
    const double d1 = (1.0 - w) * d1_0 + w * d1_1;
    const double d2 = (1.0 - w) * d2_0 + w * d2_1;
    
    return std::make_tuple(d0, d1, d2);
}


std::tuple<double, double, double> Interpolator::interpolate(const double t)
{   
    // Handle out of bounds - return boundary value with zero derivatives
    if (t <= this->t_data.front())
    {
        return std::make_tuple(this->x_data.front(), 0.0, 0.0);
    }
    else if (t >= this->t_data.back())
    {
        return std::make_tuple(this->x_data.back(), 0.0, 0.0);
    }
    
    // Use blended quartic Lagrange interpolation with smooth stencil transitions
    return lagrange_quartic_blended(t, this->t_data, this->x_data);
}*/
/*

*/



// Helper function to print expected CSV format
static void print_interpolator_format_help()
{
    std::cout << colors::bold << "Expected CSV format:" << colors::reset << "\n"
              << " * Header: 1 line (e.g., 't [s], R [m], R_dot [m/s]' or 't [s], p [Pa]')\n"
              << " * Seperator: comma (',')\n"
              << " * First column: time values in seconds (strictly monotonically increasing)\n"
              << " * Second column: corresponding data values (bubble radius in meters or ambient pressure in pascals)\n"
			  << " * Third column: corresponding data values (time derivative of bubble radius in meters per second)\n"
              << " * At least 10 data points required\n"
              << " * Example:\n"
              << "      t [s],R [m], R_dot [m]\n"
              << "      0.0,1.0e-5,0.0\n"
              << "      1.0e-6,1.001e-5,0.01\n"
              << "      2.0e-6,1.002e-5,0.01\n"
              << "      ...\n" << std::endl;
}


Interpolator::Interpolator()
{
    this->error_ID = ErrorHandler::no_error;
}


Interpolator::Interpolator(std::string csv_path)
{
    Timer timer; timer.start();
    this->error_ID = ErrorHandler::no_error;
    
    // Convert to absolute path
    std::filesystem::path abs_path = std::filesystem::absolute(csv_path);
    this->filename = abs_path.string();
    
    // Replace \ and \\ with / for consistent path representation
    for (size_t i = 0; i < this->filename.length(); ++i)
    {
        if (this->filename[i] == '\\')
        {
            this->filename[i] = '/';
        }
    }
    size_t pos = 0;
    while ((pos = this->filename.find("//", pos)) != std::string::npos)
    {
        this->filename.replace(pos, 2, "/");
    }
    
    // Check if file exists
    std::ifstream file(csv_path);
    if (!file.is_open())
    {
        this->error_ID = LOG_ERROR("Failed to open CSV file: " + csv_path);
        return;
    }
    
    std::string line;
    std::vector<double> t_temp, x_temp, v_temp;
    
    // Skip header row
    if (!std::getline(file, line))
    {
        print_interpolator_format_help();
        this->error_ID = LOG_ERROR("CSV file is empty: " + this->filename);
        file.close();
        return;
    }
    
    // Read data rows
    size_t row_count = 2;
    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string col1, col2, col3;
        
        // Read columns
        if (!std::getline(ss, col1, ',') || !std::getline(ss, col2, ',') || !std::getline(ss, col3, ',')) 
        {
            print_interpolator_format_help();
            this->error_ID = LOG_ERROR(
                "Not enough columns in row " + std::to_string(row_count) + " (\"" + line + "\") in \"" + this->filename + "\""
            );
            file.close();
            return;
        }
        
        // Parse doubles
        double t_val, x_val, v_val;
        try
        {
            t_val = std::stod(col1);
            x_val = std::stod(col2);
			v_val = std::stod(col3);
        }
        catch (const std::exception& e)
        {
            print_interpolator_format_help();
            this->error_ID = LOG_ERROR(
                "Failed to parse doubles in row " + std::to_string(row_count) + " (\"" + line + "\") in \"" + this->filename + "\""
            );
            file.close();
            return;
        }
        
        t_temp.push_back(t_val);
        x_temp.push_back(x_val);
		v_temp.push_back(v_val);
        row_count++;
    }
    file.close();
    
    // Check minimum data points
    if (t_temp.size() < 10)
    {
        print_interpolator_format_help();
        this->error_ID = LOG_ERROR(
            "At least 10 data points are required, but only got " + std::to_string(t_temp.size()) + " in \"" + this->filename + "\""
        );
        return;
    }

    // Remove exponentially growing initial steps, if present
    const double max_spacing_ratio = 50.0;
    double x_min = x_temp[0], x_max = x_temp[0];
    for (size_t i = 1; i < x_temp.size(); ++i)
    {
        if (x_temp[i] < x_min) x_min = x_temp[i];
        if (x_temp[i] > x_max) x_max = x_temp[i];
    }
    const double data_range = std::abs(x_max - x_min);
    const double quiet_tol = 1e-8 * data_range;
    
    size_t first_good = 0;
    for (size_t i = 0; i < t_temp.size() - 4; ++i)
    {
        // Calculate 4 gaps for a 5-point stencil
        double dt_min = t_temp[i + 1] - t_temp[i];
        double dt_max = dt_min;
        for (size_t j = i + 1; j < i + 4; ++j)
        {
            double dt = t_temp[j + 1] - t_temp[j];
            if (dt < dt_min) dt_min = dt;
            if (dt > dt_max) dt_max = dt;
        }
        
        const double ratio = (dt_min > 0) ? (dt_max / dt_min) : 1e10;
        const double val_change = std::abs(x_temp[i + 4] - x_temp[i]);
        
        if (ratio > max_spacing_ratio && val_change < quiet_tol)
        {
            first_good = i + 1;  // drop this point
        }
        else
        {
            break;
        }
    }
    
    if (first_good > 0)
    {
        LOG_ERROR(
            Error::severity::info,
            Error::type::preprocess,
            "Removed " + std::to_string(first_good) + " ill-conditioned startup points from \"" + this->filename + "\"" 
        );
        
        // Shift data to remove startup points
        std::vector<double> t_filtered(t_temp.begin() + first_good, t_temp.end());
        std::vector<double> x_filtered(x_temp.begin() + first_good, x_temp.end());
		std::vector<double> v_filtered(v_temp.begin() + first_good, v_temp.end());
        t_temp = std::move(t_filtered);
        x_temp = std::move(x_filtered);
		v_temp = std::move(v_filtered);
    }
    
    // Remove non-monotonous points (strict inequality to prevent duplicates)
    size_t removed_count = 0;
    this->t_data.reserve(t_temp.size());
    this->x_data.reserve(t_temp.size());
	this->v_data.reserve(t_temp.size());
    this->t_data.push_back(t_temp[0]);
    this->x_data.push_back(x_temp[0]);
	this->v_data.push_back(v_temp[0]);
    
    for (size_t i = 1; i < t_temp.size(); ++i)
    {
        if (t_temp[i] > this->t_data.back())  // Changed from >= to >
        {
            this->t_data.push_back(t_temp[i]);
            this->x_data.push_back(x_temp[i]);
			this->v_data.push_back(v_temp[i]);
        }
        else
        {
            removed_count++;
        }
    }
    
    if (removed_count > 0)
    {
        LOG_ERROR(
            Error::severity::info,
            Error::type::preprocess,
            "Removed " + std::to_string(removed_count) + " not strictly monotonious points from \"" + this->filename + "\"" 
        );
    }
    
    // Log success with timing
    const double runtime = timer.lap();
    LOG_ERROR(
        Error::severity::info,
        Error::type::preprocess,
        "Loaded " + std::to_string(this->t_data.size()) + " data points from \"" + this->filename + "\" in " + Timer::format_time(runtime)
    );
}


// Lagrange basis polynomial L_i(t) for i-th point among 5 stencil points
inline double lagrange_basis_L(const double t, const double ti, const double tj, const double tk, const double tl, const double tm)
{
    return (t - tj) * (t - tk) * (t - tl) * (t - tm) /
           ((ti - tj) * (ti - tk) * (ti - tl) * (ti - tm));
}


// First derivative of Lagrange basis polynomial L_i'(t)
inline double lagrange_basis_dL(const double t, const double ti, const double tj, const double tk, const double tl, const double tm)
{
    const double denom = (ti - tj) * (ti - tk) * (ti - tl) * (ti - tm);
    return ((t - tk) * (t - tl) * (t - tm) +
            (t - tj) * (t - tl) * (t - tm) +
            (t - tj) * (t - tk) * (t - tm) +
            (t - tj) * (t - tk) * (t - tl)) / denom;
}


// Second derivative of Lagrange basis polynomial L_i''(t)
inline double lagrange_basis_ddL(const double t, const double ti, const double tj, const double tk, const double tl, const double tm)
{
    const double denom = (ti - tj) * (ti - tk) * (ti - tl) * (ti - tm);
    return 2.0 * ((t - tl) * (t - tm) + (t - tk) * (t - tm) + (t - tk) * (t - tl) +
                  (t - tj) * (t - tm) + (t - tj) * (t - tl) + (t - tj) * (t - tk)) / denom;
}


// Quartic Lagrange polynomial evaluation at point t using 5-point stencil
// Returns interpolated value P(t), first derivative P'(t), and second derivative P''(t)
// Uses explicit formulation of Lagrange basis polynomials for clarity and performance
// Implements local time normalization to [0,1] for improved numerical conditioning
inline std::tuple<double, double, double> lagrange_quartic_raw(
    const double t,
    const std::span<const double> t_pts,
    const std::span<const double> x_pts,
	const std::span<const double> v_pts
)
{
    const double t0 = t_pts[0], t1 = t_pts[1], t2 = t_pts[2], t3 = t_pts[3], t4 = t_pts[4];
    const double y0 = x_pts[0], y1 = x_pts[1], y2 = x_pts[2], y3 = x_pts[3], y4 = x_pts[4];
	const double v0 = v_pts[0], v1 = v_pts[1], v2 = v_pts[2], v3 = v_pts[3], v4 = v_pts[4];
    
    // Local time normalization: map [t_min, t_max] → [0, 1]
    // This dramatically improves numerical conditioning when timesteps vary by orders of magnitude
    const double t_min = t0;  // First point in stencil
    const double t_max = t4;  // Last point in stencil
    const double dt_range = t_max - t_min;
    
    // Normalize query point and stencil points to [0, 1]
    const double t_norm = (t - t_min) / dt_range;
    const double t0_norm = 0.0;  // (t0 - t_min) / dt_range = 0
    const double t4_norm = 1.0;  // (t4 - t_min) / dt_range = 1
    const double t1_norm = (t1 - t_min) / dt_range;
    const double t2_norm = (t2 - t_min) / dt_range;
    const double t3_norm = (t3 - t_min) / dt_range;
    
    // Interpolated value: P(t) = Σ y_i * L_i(t)
    // Value is scale-invariant, so no adjustment needed
    const double d0  = y0 * lagrange_basis_L(t_norm, t0_norm, t1_norm, t2_norm, t3_norm, t4_norm) + 
                       y1 * lagrange_basis_L(t_norm, t1_norm, t0_norm, t2_norm, t3_norm, t4_norm) +
                       y2 * lagrange_basis_L(t_norm, t2_norm, t0_norm, t1_norm, t3_norm, t4_norm) + 
                       y3 * lagrange_basis_L(t_norm, t3_norm, t0_norm, t1_norm, t2_norm, t4_norm) +
                       y4 * lagrange_basis_L(t_norm, t4_norm, t0_norm, t1_norm, t2_norm, t3_norm);
    
    // First derivative: P'(t) = Σ y_i * L_i'(t)
    // Chain rule: dP/dt_actual = dP/dt_norm × dt_norm/dt_actual = dP/dt_norm / dt_range
    const double d1_norm = v0 * lagrange_basis_L(t_norm, t0_norm, t1_norm, t2_norm, t3_norm, t4_norm) + 
                           v1 * lagrange_basis_L(t_norm, t1_norm, t0_norm, t2_norm, t3_norm, t4_norm) +
                           v2 * lagrange_basis_L(t_norm, t2_norm, t0_norm, t1_norm, t3_norm, t4_norm) + 
                           v3 * lagrange_basis_L(t_norm, t3_norm, t0_norm, t1_norm, t2_norm, t4_norm) +
                           v4 * lagrange_basis_L(t_norm, t4_norm, t0_norm, t1_norm, t2_norm, t3_norm);
	const double d1 = d1_norm;
	/*const double d1_norm = y0 * lagrange_basis_dL(t_norm, t0_norm, t1_norm, t2_norm, t3_norm, t4_norm) + 
                           y1 * lagrange_basis_dL(t_norm, t1_norm, t0_norm, t2_norm, t3_norm, t4_norm) +
                           y2 * lagrange_basis_dL(t_norm, t2_norm, t0_norm, t1_norm, t3_norm, t4_norm) + 
                           y3 * lagrange_basis_dL(t_norm, t3_norm, t0_norm, t1_norm, t2_norm, t4_norm) +
                           y4 * lagrange_basis_dL(t_norm, t4_norm, t0_norm, t1_norm, t2_norm, t3_norm);
    const double d1 = d1_norm / dt_range;*/
    
    // Second derivative: P''(t) = Σ y_i * L_i''(t)
    // Chain rule: d²P/dt_actual² = d²P/dt_norm² × (dt_norm/dt_actual)² = d²P/dt_norm² / dt_range²
    const double d2_norm = y0 * lagrange_basis_ddL(t_norm, t0_norm, t1_norm, t2_norm, t3_norm, t4_norm) + 
                           y1 * lagrange_basis_ddL(t_norm, t1_norm, t0_norm, t2_norm, t3_norm, t4_norm) +
                           y2 * lagrange_basis_ddL(t_norm, t2_norm, t0_norm, t1_norm, t3_norm, t4_norm) + 
                           y3 * lagrange_basis_ddL(t_norm, t3_norm, t0_norm, t1_norm, t2_norm, t4_norm) +
                           y4 * lagrange_basis_ddL(t_norm, t4_norm, t0_norm, t1_norm, t2_norm, t3_norm);
    const double d2 = d2_norm / (dt_range * dt_range);
    
    return std::make_tuple(d0, d1, d2);
}


// Blended quartic Lagrange interpolation with smooth C^2 stencil transitions
// Addresses the stencil discontinuity problem: when the interpolation point moves from one
// 5-point stencil to another, abrupt changes can occur. This function smoothly blends two
// overlapping stencils using Perlin smootherstep weighting (5th-order polynomial) to ensure
// continuity in value, first derivative, and second derivative across stencil boundaries.
// The transition occurs in the overlap region between adjacent stencils.
std::tuple<double, double, double> lagrange_quartic_blended(
    const double t,
    const std::vector<double>& t_data,
    const std::vector<double>& x_data,
	const std::vector<double>& v_data
)
{
    const size_t n = t_data.size();
    
    // Binary search to locate interval: t[j] <= t < t[j+1]
    size_t low = 0, high = n - 1;
    while (low <= high)
    {
        size_t mid = low + (high - low) / 2;
        if (t_data[mid] <= t)
        {
            low = mid + 1;
        }
        else
        {
            high = mid - 1;
        }
    }
    size_t j = high;
    j = std::min(std::max(j, static_cast<size_t>(0)), n - 2);
    
    // Primary stencil i0, centered around j
    size_t i0 = std::min(std::max((int)j - 2, 0), (int)(n - 5));
    
    // Shifted stencil i0b (one step to the right)
    size_t i0b = std::min(i0 + 1, n - 5);
    
    // EDGE CASE 1: Stencils are identical (happens near right boundary)
    // When i0 >= n-4, there's no room for a second stencil: i0b = min(i0+1, n-5) = n-5 = i0
    // Example: n=10 → max stencil start is n-5=5. If i0=5, then i0b=min(6,5)=5 → same stencil
    // This happens when we're in the last 5 points of data. No blending possible, use single stencil.
    if (i0 == i0b)
    {
        std::span<const double> t_pts(t_data.data() + i0, 5);
        std::span<const double> x_pts(x_data.data() + i0, 5);
		std::span<const double> v_pts(v_data.data() + i0, 5);
        return lagrange_quartic_raw(t, t_pts, x_pts, v_pts);
    }
    
    // Define blending transition interval [t_lo, t_hi]
    // Stencil A: [i0, i0+1, i0+2, i0+3, i0+4]
    // Stencil B: [i0b, i0b+1, i0b+2, i0b+3, i0b+4] = [i0+1, i0+2, i0+3, i0+4, i0+5]
    // Overlap points: i0+2, i0+3, i0+4 (shared by both stencils)
    const double t_lo = t_data[i0 + 2];  // Start of transition zone
    const double t_hi = t_data[i0 + 3];  // End of transition zone
    
    // Normalize position within transition zone to s ∈ [0, 1]
    double s = (t - t_lo) / (t_hi - t_lo);
    s = std::min(std::max(s, 0.0), 1.0);
    
    // EDGE CASE 2: Query point at or before transition start (s = 0)
    // When t ≤ t_lo, we're entirely in stencil A's domain before the blend zone
    // Blending weight w would be 0 anyway (since w(0)=0), so skip computation and use stencil A
    if (s <= 0.0)
    {
        std::span<const double> t_pts(t_data.data() + i0, 5);
        std::span<const double> x_pts(x_data.data() + i0, 5);
		std::span<const double> v_pts(v_data.data() + i0, 5);
        return lagrange_quartic_raw(t, t_pts, x_pts, v_pts);
    }
    
    // EDGE CASE 3: Query point at or after transition end (s = 1)
    // When t ≥ t_hi, we're entirely in stencil B's domain after the blend zone
    // Blending weight w would be 1 anyway (since w(1)=1), so skip computation and use stencil B
    if (s >= 1.0)
    {
        std::span<const double> t_pts(t_data.data() + i0b, 5);
        std::span<const double> x_pts(x_data.data() + i0b, 5);
		std::span<const double> v_pts(v_data.data() + i0b, 5);
        return lagrange_quartic_raw(t, t_pts, x_pts, v_pts);
    }
    
    // Perlin smootherstep: C^2 smooth weight function
    // w = s^3(6s^2 - 15s + 10)
    const double w = s * s * s * (s * (s * 6.0 - 15.0) + 10.0);
    
    // Evaluate both stencils
    std::span<const double> t_pts0(t_data.data() + i0, 5);
    std::span<const double> x_pts0(x_data.data() + i0, 5);
	std::span<const double> v_pts0(v_data.data() + i0, 5);
    auto [d0_0, d1_0, d2_0] = lagrange_quartic_raw(t, t_pts0, x_pts0, v_pts0);
    
    std::span<const double> t_pts1(t_data.data() + i0b, 5);
    std::span<const double> x_pts1(x_data.data() + i0b, 5);
	std::span<const double> v_pts1(v_data.data() + i0b, 5);
    auto [d0_1, d1_1, d2_1] = lagrange_quartic_raw(t, t_pts1, x_pts1, v_pts1);
    
    // Blend results
    const double d0 = (1.0 - w) * d0_0 + w * d0_1;
    const double d1 = (1.0 - w) * d1_0 + w * d1_1;
    const double d2 = (1.0 - w) * d2_0 + w * d2_1;
    
    return std::make_tuple(d0, d1, d2);
}


std::tuple<double, double, double> Interpolator::interpolate(const double t)
{   
    // Handle out of bounds - return boundary value with zero derivatives
    if (t <= this->t_data.front())
    {
        return std::make_tuple(this->x_data.front(), 0.0, 0.0);
    }
    else if (t >= this->t_data.back())
    {
        return std::make_tuple(this->x_data.back(), 0.0, 0.0);
    }
    
    // Use blended quartic Lagrange interpolation with smooth stencil transitions
    return lagrange_quartic_blended(t, this->t_data, this->x_data, this->v_data);
}
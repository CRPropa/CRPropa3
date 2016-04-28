#ifndef CRPROPA_COMMON_H
#define CRPROPA_COMMON_H

#include <string>
#include <vector>

/**
 @file
 @brief Common helper functions
 */

namespace crpropa {

// Returns the full path to a CRPropa data file
std::string getDataPath(std::string filename);

// Returns a certain digit from a given integer
inline int digit(const int& value, const int& d) {
	return (value % (d * 10)) / d;
}

// Return value xclip which is the closest to x, so that lower <= xclip <= upper
template <typename T>
T clip(const T& x, const T& lower, const T& upper) {
	return std::max(lower, std::min(x, upper));
}

// Perform linear interpolation on a set of n tabulated data points X[0 .. n-1] -> Y[0 .. n-1]
// Returns Y[0] if x < X[0] and Y[n-1] if x > X[n-1]
double interpolate(double x, const std::vector<double>& X,
		const std::vector<double>& Y);


// Perform bilinear interpolation on a set of (n,m) tabulated data points X[0 .. n-1], Y[0 .. m-1] -> Z[0.. n-1*m-1]
// Returns 0 if x < X[0] or x > X[n-1] or y < Y[0] or y > Y[m-1]
double interpolate2d(double x, double y, const std::vector<double>& X,
		const std::vector<double>& Y, const std::vector<double>& Z);

// Perform linear interpolation on equidistant tabulated data
// Returns Y[0] if x < lo and Y[n-1] if x > hi
double interpolateEquidistant(double x, double lo, double hi,
		const std::vector<double>& Y);

} // namespace crpropa

#endif // CRPROPA_COMMON_H

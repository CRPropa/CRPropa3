#ifndef CRPROPA_COMMON_H
#define CRPROPA_COMMON_H

#include <string>
#include <vector>

/**
 @file
 @brief Common helper functions
 */

namespace crpropa {
/**
 * \addtogroup Core
 * @{
 */

// Returns the full path to a CRPropa data file
std::string getDataPath(std::string filename);

// Returns the install prefix
std::string getInstallPrefix();

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

// Find index of value in a sorted vector X that is closest to x
size_t closestIndex(double x, const std::vector<double> &X);
/** @}*/


// pow implementation as template for integer exponents pow_integer<2>(x)
// evaluates to x*x
template <unsigned int exponent>
inline double pow_integer(double base)
{
  return pow_integer<(exponent >> 1)>(base*base) * (((exponent & 1) > 0) ? base : 1);
}

template <>
inline double pow_integer<0>(double base)
{
  return 1;
}



} // namespace crpropa

#endif // CRPROPA_COMMON_H

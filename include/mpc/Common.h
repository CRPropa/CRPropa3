#ifndef COMMON_H_
#define COMMON_H_

#include <string>
#include <vector>

namespace mpc {

// Photon fields
enum PhotonField {
	CMB, IRB, CMB_IRB
};

// Returns overall photon field scaling factor at redshift z
double photonFieldScaling(int photonField, double z);

// Returns the full path to a mpc data file
std::string getDataPath(std::string filename);

// Returns a certain digit from a given integer
inline int digit(const int& value, const int& d) {
	return (value % (d * 10)) / d;
}

// Perform linear interpolation on a set of n tabulated data points X[0 .. n-1] -> Y[0 .. n-1]
// Returns Y[0] if x < X[0] and Y[n-1] if x > X[n-1]
double interpolate(double x, const std::vector<double>& X,
		const std::vector<double>& Y);

// Perform linear interpolation on equidistant tabulated data
// Returns Y[0] if x < lo and Y[n-1] if x > hi
double interpolateEquidistant(double x, double lo, double hi,
		const std::vector<double>& Y);

} // namespace mpc

#endif /* COMMON_H_ */

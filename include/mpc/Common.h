#ifndef COMMON_H_
#define COMMON_H_

#include <string>
#include <vector>

namespace mpc {

// Photon fields
enum PhotonField {
	CMB, IRB, CMB_IRB
};

// Returns the full path to a mpc data file
std::string getDataPath(std::string filename);

// Returns a certain digit from a given integer
inline int digit(const int &value, const int &d) {
	return (value % (d * 10)) / d;
}

// Perform linear interpolation
double interpolate(double x, const std::vector<double> &X,
		const std::vector<double> &Y);

// Perform linear interpolation on equidistant tabulated data
double interpolateEquidistant(double x, double lo, double hi,
		const std::vector<double> &Y);

} // namespace mpc

#endif /* COMMON_H_ */

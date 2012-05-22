#ifndef COMMON_H_
#define COMMON_H_

#include <string>

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
double interpolate(const double x, const double *xD, const double *yD);

// Perform linear interpolation on equidistant tabulated data
double interpolateEquidistant(const double x, const double xLo, const double dx,
		const double *yD);

} // namespace mpc

#endif /* COMMON_H_ */

#ifndef COMMON_H_
#define COMMON_H_

#include <string>

namespace mpc {

// Returns the full path to a mpc data file
std::string getDataPath(std::string filename);

// Returns a certain digit from a given integer
inline int digit(const int &value, const int &d) {
	return (value % (d * 10)) / d;
}

// Photon fields
enum PhotonField {
	CMB, IR, CMBIR
};

} // namespace mpc

#endif /* COMMON_H_ */

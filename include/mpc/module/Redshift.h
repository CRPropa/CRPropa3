#ifndef REDSHIFT_H_
#define REDSHIFT_H_

#include "mpc/Module.h"

#include <vector>

namespace mpc {

/**
 @class SimpleRedshift
 @brief Redshift and adiabatic energy loss using Hubble's law and the distance to an observer.

 Using Hubble's law v = H0*D and the the small redshift approximation v ~ c*z, the redshift is calculated as z ~ H0*D/c.
 It is assumed that the particle starts with a redshift z, that corresponds to its distance to the observer.
 This redshift is reduced with shrinking distance to the observer, and the particle loses energy accordingly.
 Redshift and particle energy are not changed if the distance to the observer grows.
 */
class SimpleRedshift : public Module {
private:
	Vector3d observer; // observer position (z = 0)
	double h; // dimension-free Hubble constant, H0 = h * 100 km/s/Mpc

public:
	SimpleRedshift(Vector3d observer, double h = 0.7);
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

///**
// @class Redshift
// @brief Redshift and adiabatic energy loss
// */
//class Redshift {
//private:
//	std::vector<double> z;
//	std::vector<double> D;
//
//public:
////	double h; // dimension-free Hubble constant, H0 = h * 100 km/s/Mpc
////	double omegaM; // density parameter
////	double omegaL; // vacuum energy parameter
//	Redshift(double h = 0.7, double omegaM = 0.3, double omegaL = 0.7);
//	void process(Candidate *candidate) const;
//	std::string getDescription() const;
//};

} // namespace mpc

#endif /* REDSHIFT_H_ */

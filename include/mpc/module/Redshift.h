#ifndef MPC_REDSHIFT_H_
#define MPC_REDSHIFT_H_

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

/**
 @class Redshift
 @brief Calculation of cosmological redshift and adiabatic energy loss
 */
class Redshift : public Module {
private:
	double H0; // Hubble rate at z=0 in [1/s], H0 = h * 100 km/s/Mpc
	double omegaM; // density parameter
	double omegaL; // vacuum energy parameter

	static const int n = 1000;
	static const double zmin = 0.0001;
	static const double zmax = 100;

	std::vector<double> Z; // redshift
	std::vector<double> D; // comoving distance [m]

public:
	Redshift(double h = 0.7, double omegaM = 0.3, double omegaL = 0.7);
	void init();
	double getRedshift(double distance) const;
	double getHubbleRate(double redshift) const;
	double getDistance(double redshift) const;
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namespace mpc

#endif /* MPC_REDSHIFT_H_ */

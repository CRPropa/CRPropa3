#include "mpc/module/Redshift.h"
#include "mpc/Units.h"
#include "mpc/Common.h"

#include <sstream>
#include <stdexcept>

namespace mpc {

SimpleRedshift::SimpleRedshift(Vector3d observer, double h) :
		observer(observer), h(h) {
}

void SimpleRedshift::process(Candidate *c) const {
	double d = c->current.getPosition().getDistanceTo(observer);
	double z = h * 1e5 / Mpc * d / c_light;
	double dz = c->getRedshift() - z;

	if (dz < 0)
		return; // do nothing if the new redshift is higher

	c->setRedshift(z);
	double E = c->current.getEnergy();
	c->current.setEnergy(E * (1 - dz / (1 + z))); // using dE/dz = E/(1+z)
}

std::string SimpleRedshift::getDescription() const {
	std::stringstream ss;
	ss << "Simple redshift: observer " << observer / Mpc << "Mpc, h = " << h;
	return ss.str();
}

Redshift::Redshift(double h, double m, double l) :
		H0(h * 1e5 / Mpc), omegaM(m), omegaL(l) {
	Z.resize(n);
	D.resize(n);
	std::vector<double> H(n);

	Z[0] = 0;
	D[0] = 0;
	H[0] = H0;

	// Relation between comoving distance r and redshift z (cf. J.A. Peacock, Cosmological physics, p. 89 eq. 3.76)
	// dr = c / H(z) dz, integration using midpoint rule
	double dlz = log10(zmax) - log10(zmin);
	double dz;
	for (int i = 1; i < n; i++) {
		Z[i] = zmin * pow(10, i * dlz / (n - 1)); // logarithmic even spacing
		dz = (Z[i] - Z[i - 1]);
		H[i] = hubbleRate(Z[i]);
		D[i] = D[i - 1]	+ c_light * (1 / H[i] + 1 / H[i - 1]) / 2 * dz;
	}
}

double Redshift::hubbleRate(double z) const {
	return H0 * sqrt(omegaL + omegaM * pow(1 + z, 3));
}

double Redshift::comovingDistance2Redshift(double d) const {
	if (d < 0)
		throw std::runtime_error("Redshift: d < 0");
	if (d > D[n - 1])
		throw std::runtime_error("Redshift: d > dmax");
	return interpolate(d, D, Z);
}

double Redshift::redshift2ComovingDistance(double z) const {
	if (z < 0)
		throw std::runtime_error("Redshift: z < 0");
	if (z > zmax)
		throw std::runtime_error("Redshift: z > zmax");
	double d = interpolate(z, Z, D);
	return d;
}

void Redshift::process(Candidate *c) const {
	double z = c->getRedshift();
	if (z == 0)
		return; // nothing to do, redshift can't get smaller

	// use small redshift approximation:  dz = H(z) / c * dr
	double dz = hubbleRate(z) / c_light * c->getCurrentStep();

	// update redshift
	c->setRedshift(z - dz);

	// adiabatic energy loss: dE / dz = E/(1+z)
	double E = c->current.getEnergy();
	c->current.setEnergy(E * (1 - dz / (1 + z)));
}

std::string Redshift::getDescription() const {
	std::stringstream ss;
	ss << "Redshift: " << "H0 = " << H0 / 1e3 * Mpc << " km/s/Mpc, ";
	ss << "OmegaL = " << omegaL << ", OmegaM = " << omegaM;
	return ss.str();
}

} // namespace mpc

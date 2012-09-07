#include "mpc/module/Redshift.h"
#include "mpc/Units.h"
#include "mpc/Common.h"

#include <sstream>

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
	init();
}

void Redshift::init() {
	Z.resize(n);
	D.resize(n);
	std::vector<double> H(n);

	Z[0] = 0;
	D[0] = 0;
	H[0] = H0;

	// Relation between comoving distance and redshift (see J.A. Peacock, Cosmological physics, p. 89 eq. 3.76)
	// R0 dr = c / H(z) dz
	// Integration with midpoint rule: R0 dr = (c/H(z+dz) + c/H(z)) / 2 * dz
	double dlz = log10(zmax) - log10(zmin);
	double dz, A;
	for (int i = 1; i < n; i++) {
		Z[i] = zmin * pow(10, i * dlz / (n - 1)); // logarithmic even spacing
		H[i] = getHubbleRate(Z[i]);
		D[i] = D[i - 1]
				+ (Z[i] - Z[i - 1]) * c_light * (1 / H[i] + 1 / H[i - 1]) / 2;
	}
}

double Redshift::getHubbleRate(double z) const {
	return H0 * sqrt(omegaL + omegaM * pow(1 + z, 3));
}

double Redshift::getRedshift(double d) const {
	if (d < 0)
		throw std::runtime_error("Redshift: d < 0");
	if (d > D[n - 1])
		throw std::runtime_error("Redshift: d > dmax");
	return interpolate(d, D, Z);
}

double Redshift::getDistance(double z) const {
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

	// use small step approximation to calculate redshift change
	// dz = H(z) / c * dr
	// dE/dz = E/(1+z)
	double dz = getHubbleRate(z) / c_light * c->getCurrentStep();

	// update candidate
	c->setRedshift(z - dz);
	double E = c->current.getEnergy();
	c->current.setEnergy(E * (1 - dz / (1 + z)));
}

std::string Redshift::getDescription() const {
	std::stringstream ss;
	ss << "Redshift for cosmological parameters: ";
	ss << "H0 = " << H0 * 1e3 / Mpc << " km/s/Mpc, ";
	ss << "omegaL = " << omegaL << ", omegaM = " << omegaM;
	return ss.str();
}

} // namespace mpc

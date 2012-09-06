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
	double dz = z - c->getRedshift();

	if (dz > 0)
		return; // do nothing if the new redshift is higher

	c->setRedshift(z);

	double E = c->current.getEnergy();
	c->current.setEnergy(E * (1 + dz) / (1 + z)); // using dE/dz=E/(1+z)
}

std::string SimpleRedshift::getDescription() const {
	std::stringstream ss;
	ss << "Simple redshift: observer " << observer / Mpc << "Mpc, h = " << h;
	return ss.str();
}

Redshift::Redshift(double h, double m, double l) {
	H0 = h * 1e5 / Mpc;
	omegaM = m;
	omegaL = l;

	std::stringstream ss;
	ss << "Redshift for cosmological parameters: h = " << h << ", omegaL = "
			<< omegaL << ", omegaM = " << omegaM;
	setDescription(ss.str());
}

void Redshift::process(Candidate *c) const {
	double z = c->getRedshift();
	if (z == 0)
		return; // nothing to do

	// small step approximation: dz = H(z) / c * dr
	// an analytical form for the integral would be better
	double H = H0 * sqrt(omegaL + omegaM * pow(1 + z, 3));
	double dz = H / c_light * c->getCurrentStep();
	dz = std::min(z, dz); // dz cannot be larger than z

	c->setRedshift(z - dz);
	double E = c->current.getEnergy();
	c->current.setEnergy(E * (1 - dz) / (1 + z)); // using dE/dz=E/(1+z)
}

} // namespace mpc

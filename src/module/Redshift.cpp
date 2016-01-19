#include "crpropa/module/Redshift.h"
#include "crpropa/Units.h"
#include "crpropa/Cosmology.h"

#include <limits>

namespace crpropa {

void Redshift::process(Candidate *c) const {
	double z = c->getRedshift();

	// check if z = 0
	if (z <= std::numeric_limits<double>::min())
		return;

	// use small step approximation:  dz = H(z) / c * ds
	double dz = hubbleRate(z) / c_light * c->getCurrentStep();

	// prevent dz > z
	dz = std::min(dz, z);

	// update redshift
	c->setRedshift(z - dz);

	// adiabatic energy loss: dE / dz = E / (1 + z)
	double E = c->current.getEnergy();
	c->current.setEnergy(E * (1 - dz / (1 + z)));
}

std::string Redshift::getDescription() const {
	std::stringstream s;
	s << "Redshift: h0 = " << hubbleRate() / 1e5 * Mpc << ", omegaL = "
			<< omegaL() << ", omegaM = " << omegaM();
	return s.str();
}

void FutureRedshift::process(Candidate *c) const {
	double z = c->getRedshift();

	// check if z = -1
	if (z <= -1)
		return;

	// use small step approximation:  dz = H(z) / c * ds
	double dz = hubbleRate(z) / c_light * c->getCurrentStep();

	// update redshift
	c->setRedshift(z - dz);

	// adiabatic energy loss: dE / dz = E / (1 + z)
	double E = c->current.getEnergy();
	c->current.setEnergy(E * (1 - dz / (1 + z)));
}

std::string FutureRedshift::getDescription() const {
	std::stringstream s;
	s << "FutureRedshift: h0 = " << hubbleRate() / 1e5 * Mpc << ", omegaL = "
			<< omegaL() << ", omegaM = " << omegaM();
	return s.str();
}

} // namespace crpropa

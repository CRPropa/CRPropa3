#include "mpc/module/Redshift.h"
#include "mpc/Cosmology.h"

namespace mpc {

Redshift::Redshift() {
	setDescription("Redshift");
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

} // namespace mpc

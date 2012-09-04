#include "mpc/module/Redshift.h"
#include "mpc/Units.h"

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

//Redshift::Redshift(double h, double omegaM, double omegaL) {
//	double H0 = h * 1e5 / Mpc;
//
//	const int n = 1000;
//	double zMin = 0.0001;
//	double zMax = 100;
//
//	std::vector<double> H;
//
//	z.resize(n);
//	H.resize(n);
//	D.resize(n);
//
//	z[0] = 0;
//	H[0] = H0;
//	D[0] = 0;
//
//	// Relation between comoving distance and redshift (see J.A. Peacock, Cosmological physics, p. 89 eq. 3.76)
//	// R0 dr = c / H(z) dz
//	// H(z) = H0 sqrt(omegaL + omegaM (1 + z)^3)
//	// Integration with midpoint rule.
//	for (int i = 1; i < n; i++) {
//		z[i] = pow(10, zMin + (zMax - zMin) * i / (n - 1));
//		H[i] = H0 * sqrt(omegaL + omegaM * pow(1 + z[i], 3));
//		D[i] = D[i - 1]
//				+ c_light * (z[i] - z[i - 1]) * (1 / H[i - 1] + 1 / H[i]) / 2;
//	}
//}
//
//void Redshift::process(Candidate *c) const {
//}
//
//std::string Redshift::getDescription() const {
//	std::stringstream ss;
//	ss << "Redshift: h = ";
//	return ss.str();
//}

//void TParticlePropa::RedshiftLoss() {
//	if (_fpUniverse->Type() == UNIVERSE_ENV1D) {
//		// Using the relations dE/dz=E/(1+z) and dz/dt=H_0(1+z)E(z)
//		double lLoss = _fpUniverse->H0() * _fTimeStep * sqrt(_fpUniverse->OmegaLambda() + _fpUniverse->OmegaM() * pow(1 + _fRedshift, 3));
//		_fEnergy *= (1 - lLoss);
//		_fMomentum.setMag(_fEnergy * inv_c_light);
//	}
//}

} // namespace mpc

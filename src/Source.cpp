#include "mpc/Source.h"

#include "mpc/Random.h"

namespace mpc {

BasicSource::BasicSource() :
		position(Vector3d(0, 0, 0)), id(0), index1(-1), index2(-1), breakpoint(
				10 * EeV), Emin(5 * EeV), Emax(1000 * EeV) {

}

BasicSource::BasicSource(const Vector3d &sposition, int type, double Emin,
		double Emax, double index1, double index2, double breakpoint) :
		position(position), id(id), index1(index1), index2(index2), breakpoint(
				breakpoint), Emin(Emin), Emax(Emax) {
}

void BasicSource::prepare(ParticleState &state) const {

	state.setPosition(position);
	state.setId(id);
	double E = Random::instance().randBrokenPowerLaw(index1, index2, breakpoint,
			Emin, Emax);
	state.setEnergy(E);
	state.setDirection(Random::instance().randUnitVectorOnSphere());
}

} // namespace mpc

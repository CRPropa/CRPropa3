#include "mpc/Source.h"

#include "mpc/Random.h"

namespace mpc {

BasicSource::BasicSource(const Vector3 &sposition, int type, double Emin,
		double Emax, double index1, double index2, double breakpoint) :
		position(position), id(id), index1(index1), index2(index2), breakpoint(
				breakpoint), Emin(Emin), Emax(Emax) {
}

void BasicSource::prepare(ParticleState &state) const {
	state.setPosition(position);
	state.setId(id);
	state.setEnergy(
			Random::instance().randBrokenPowerLaw(index1, index2, breakpoint,
					Emin, Emax));
	state.setDirection(Random::instance().randUnitVectorOnSphere());
}

} // namespace mpc

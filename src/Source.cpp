#include "mpc/Source.h"
#include "mpc/Random.h"

#include "HepPID/ParticleIDMethods.hh"

#include <stdexcept>

namespace mpc {

BasicSource::BasicSource() :
		position(Vector3d(0, 0, 0)), id(0), index1(-1), index2(-1), breakpoint(
				10 * EeV), Emin(5 * EeV), Emax(1000 * EeV) {
}

BasicSource::BasicSource(const Vector3d &position, int id, double Emin,
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

CompositeSource::CompositeSource() :
		position(Vector3d(0, 0, 0)), Emin(5), Emax(1000), index(-1) {
}

CompositeSource::CompositeSource(const Vector3d& position, double Emin,
		double Emax, double index) :
		position(position), Emin(Emin), Emax(Emax), index(index) {
}

CompositeSource::Isotope::Isotope() :
		id(getNucleusId(1, 1)), abundance(1), probability(0) {
}

CompositeSource::Isotope::Isotope(int id, double abundance) :
		id(id), abundance(abundance), probability(0) {
}

void CompositeSource::addToComposition(int id, double abundance) {
	Isotope iso(id, abundance);
	composition.push_back(iso);
	normalize();
}

void CompositeSource::normalize() {
	double pSum = 0;
	double a = 1 + index;
	for (int i = 0; i < composition.size(); i++) {
		Isotope &iso = composition[i];
		int Z = HepPID::Z(iso.id);
		if (std::abs(a) < std::numeric_limits<double>::min()) // index = -1
			pSum += iso.abundance * log(Z * Emax / Emin);
		else
			// index <> -1
			pSum += iso.abundance / a * (pow(Z * Emax, a) - pow(Emin, a));
		iso.probability = pSum;
	}
	for (int i = 0; i < composition.size(); i++) {
		composition[i].probability /= pSum;
	}
}

void CompositeSource::prepare(ParticleState& state) const {
	if (composition.size() == 0)
		throw std::runtime_error("mpc::CompositeSource: No source isotope set");
	double r = Random().rand();
	int i = -1;
	while (r < composition[i + 1].probability)
		i++;
	state.setId(composition[i].id);
	double E = Random::instance().randPowerLaw(index, Emin, Emax);
	state.setEnergy(E);
	state.setPosition(position);
	state.setDirection(Random::instance().randUnitVectorOnSphere());
}

} // namespace mpc

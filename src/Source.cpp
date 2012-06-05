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

void CompositeSource::addToComposition(int id, double a) {
	isotope.push_back(id);
	abundance.push_back(a);
	probability.push_back(0);
	normalize();
}

void CompositeSource::normalize() {
	double pSum = 0;
	double a = 1 + index;
	for (int i = 0; i < isotope.size(); i++) {
		int Z = HepPID::Z(isotope[i]);
		if (std::abs(a) < std::numeric_limits<double>::min())
			pSum += abundance[i] * log(Z * Emax / Emin);
		else
			pSum += abundance[i] / a * (pow(Z * Emax, a) - pow(Emin, a));
		probability[i] = pSum;
	}
	for (int i = 0; i < probability.size(); i++) {
		probability[i] /= pSum;
	}
}

void CompositeSource::prepare(ParticleState& state) const {
	if (isotope.size() == 0)
		throw std::runtime_error("mpc::CompositeSource: No source isotope set");
	double r = Random().rand();
	int i = 0;
	while ((r > probability[i]) and (i < probability.size()))
		i++;
	state.setId(isotope[i]);
	int Z = HepPID::Z(isotope[i]);
	state.setEnergy(Random::instance().randPowerLaw(index, Emin, Z * Emax));
	state.setPosition(position);
	state.setDirection(Random::instance().randUnitVectorOnSphere());
}

} // namespace mpc

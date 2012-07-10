#include "mpc/Source.h"
#include "mpc/Random.h"

#include "HepPID/ParticleIDMethods.hh"
#include <stdexcept>

namespace mpc {

void Source::addProperty(SourceProperty* property) {
	properties.push_back(property);
}

void Source::prepare(ParticleState& state) const {
	for (int i = 0; i < properties.size(); i++)
		(*properties[i]).prepare(state);
}

SourceParticleType::SourceParticleType(int id) :
		id(id) {
}

void SourceParticleType::prepare(ParticleState &particle) const {
	particle.setId(id);
}

SourcePowerLawSpectrum::SourcePowerLawSpectrum(double Emin, double Emax,
		double index) :
		Emin(Emin), Emax(Emax), index(index) {
}

void SourcePowerLawSpectrum::prepare(ParticleState &particle) const {
	double E = Random().randPowerLaw(index, Emin, Emax);
	particle.setEnergy(E);
}

SourceComposition::SourceComposition(double Emin, double Rmax, double index) :
		Emin(Emin), Rmax(Rmax), index(index) {
}

double SourceComposition::getSpectrumIntegral(int Z) const {
	double a = 1 + index;
	double Emax = Z * Rmax;
	if (std::abs(a) < std::numeric_limits<double>::min())
		return log(Emax / Emin);
	else
		return (pow(Emax, a) - pow(Emin, a)) / a;
}

void SourceComposition::add(int id, double a) {
	isotope.push_back(id);
	int A = getMassNumberFromNucleusId(id);
	double weightedAbundance = a * pow(A, -index - 1);
	abundance.push_back(weightedAbundance);
	probability.push_back(0);
	normalize();
}

void SourceComposition::add(int A, int Z, double a) {
	add(getNucleusId(A, Z), a);
}

void SourceComposition::normalize() {
	double pSum = 0;
	for (int i = 0; i < isotope.size(); i++) {
		int Z = HepPID::Z(isotope[i]);
		pSum += abundance[i] * getSpectrumIntegral(Z);
		probability[i] = pSum;
	}
	for (int i = 0; i < probability.size(); i++) {
		probability[i] /= pSum;
	}
}

void SourceComposition::prepare(ParticleState& particle) const {
	if (isotope.size() == 0)
		throw std::runtime_error("PowerLawComposition: No source isotope set");
	double r = Random().rand();
	int i = 0;
	while ((r > probability[i]) and (i < probability.size()))
		i++;
	int id = isotope[i];
	double E = Random().randPowerLaw(index, Emin, HepPID::Z(id) * Rmax);
	particle.setId(id);
	particle.setEnergy(E);
}

SourcePosition::SourcePosition(Vector3d position) :
		position(position) {
}

void SourcePosition::prepare(ParticleState& state) const {
	state.setPosition(position);
}

SourceSphericalVolume::SourceSphericalVolume(Vector3d center, double radius) :
		center(center), radius(radius) {
}

void SourceSphericalVolume::prepare(ParticleState& particle) const {
	double r = pow(Random().rand(), 1. / 3.) * radius;
	Vector3d pos = Random().randUnitVectorOnSphere() * r;
	particle.setPosition(pos);
}

void SourceIsotropicEmission::prepare(ParticleState& particle) const {
	Vector3d dir = Random().randUnitVectorOnSphere();
	particle.setDirection(dir);
}

SourceDirection::SourceDirection(Vector3d direction) :
		direction(direction) {
}

void SourceDirection::prepare(ParticleState& particle) const {
	particle.setDirection(direction);
}

SourceEmissionCone::SourceEmissionCone(Vector3d direction, double aperture) :
		direction(direction), aperture(aperture) {
}

void SourceEmissionCone::prepare(ParticleState& particle) const {
	Vector3d dir = Random().randConeVector(direction, aperture);
	particle.setDirection(dir);
}

} // namespace mpc

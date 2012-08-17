#include "mpc/Source.h"
#include "mpc/Random.h"

#include "HepPID/ParticleIDMethods.hh"
#include <stdexcept>

namespace mpc {

void Source::addProperty(SourceProperty *property) {
	properties.push_back(property);
}

void Source::prepare(ParticleState &particle) const {
	for (int i = 0; i < properties.size(); i++)
		(*properties[i]).prepare(particle);
}

void SourceList::addSource(Source *source, double lumi) {
	sources.push_back(source);
	if (luminosities.size() > 0)
		lumi += luminosities.back();
	luminosities.push_back(lumi);
}

void SourceList::prepare(ParticleState &particle) const {
	Random &random = Random::instance();
	if (sources.size() == 0)
		throw std::runtime_error("SourceList: no sources set");
	double r = random.rand() * luminosities.back();
	int i = 0;
	while ((r > luminosities[i]) and (i < luminosities.size()))
		i++;
	(sources[i])->prepare(particle);
}

SourceParticleType::SourceParticleType(int id) :
		id(id) {
}

void SourceParticleType::prepare(ParticleState &particle) const {
	particle.setId(id);
}

SourceEnergy::SourceEnergy(double energy) :
		E(energy) {
}

void SourceEnergy::prepare(ParticleState &particle) const {
	particle.setEnergy(E);
}

SourcePowerLawSpectrum::SourcePowerLawSpectrum(double Emin, double Emax,
		double index) :
		Emin(Emin), Emax(Emax), index(index) {
}

void SourcePowerLawSpectrum::prepare(ParticleState &particle) const {
	Random &random = Random::instance();
	double E = random.randPowerLaw(index, Emin, Emax);
	particle.setEnergy(E);
}

void SourceNuclei::add(int id, double a) {
	ids.push_back(id);
	if (abundances.size() > 0)
		a += abundances.back();
	abundances.push_back(a);
}

void SourceNuclei::prepare(ParticleState &particle) const {
	if (ids.size() == 0)
		throw std::runtime_error("SourceNuclei: no nuclei set");
	Random &random = Random::instance();
	double r = random.rand() * abundances.back();
	int i = 0;
	while ((r > abundances[i]) and (i < abundances.size()))
		i++;
	particle.setId(ids[i]);
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
	Random &random = Random::instance();
	double r = random.rand();
	int i = 0;
	while ((r > probability[i]) and (i < probability.size()))
		i++;
	int id = isotope[i];
	particle.setId(id);
	particle.setEnergy(random.randPowerLaw(index, Emin, HepPID::Z(id) * Rmax));
}

SourcePosition::SourcePosition(Vector3d position) :
		position(position) {
}

void SourcePosition::prepare(ParticleState& state) const {
	state.setPosition(position);
}

void SourceMultiplePositions::add(Vector3d pos, double lumi) {
	positions.push_back(pos);
	if (luminosities.size() > 0)
		lumi += luminosities.back();
	luminosities.push_back(lumi);
}

void SourceMultiplePositions::prepare(ParticleState &particle) const {
	if (positions.size() == 0)
		throw std::runtime_error("SourceMultiplePositions: no position set");
	double r = Random().rand() * luminosities.back();
	int i = 0;
	while ((r > luminosities[i]) and (i < luminosities.size()))
		i++;
	particle.setPosition(positions[i]);
}

SourceHomogeneousSphere::SourceHomogeneousSphere(Vector3d c, double r) :
		center(c), radius(r) {
}

void SourceHomogeneousSphere::prepare(ParticleState &particle) const {
	Random &random = Random::instance();
	double r = pow(random.rand(), 1. / 3.) * radius;
	particle.setPosition(random.randUnitVectorOnSphere() * r);
}

SourceHomogeneousBox::SourceHomogeneousBox(Vector3d o, Vector3d s) :
		origin(o), size(s) {
}

void SourceHomogeneousBox::prepare(ParticleState &particle) const {
	Random &random = Random::instance();
	Vector3d pos(random.rand(), random.rand(), random.rand());
	particle.setPosition(pos * size + origin);
}

void SourceIsotropicEmission::prepare(ParticleState &particle) const {
	Random &random = Random::instance();
	particle.setDirection(random.randUnitVectorOnSphere());
}

SourceDirection::SourceDirection(Vector3d direction) :
		direction(direction) {
}

void SourceDirection::prepare(ParticleState &particle) const {
	particle.setDirection(direction);
}

SourceEmissionCone::SourceEmissionCone(Vector3d direction, double aperture) :
		direction(direction), aperture(aperture) {
}

void SourceEmissionCone::prepare(ParticleState &particle) const {
	Random &random = Random::instance();
	particle.setDirection(random.randConeVector(direction, aperture));
}

} // namespace mpc

#include "crpropa/Source.h"
#include "crpropa/Random.h"
#include "crpropa/Cosmology.h"

#include "HepPID/ParticleIDMethods.hh"
#include <stdexcept>
#include <algorithm>

namespace crpropa {

void SourceProperty::prepare(ParticleState& particle) const {
}

void SourceProperty::prepare(Candidate& candidate) const {
	ParticleState &source = candidate.source;
	prepare(source);
	candidate.created = source;
	candidate.current = source;
	candidate.previous = source;
}

void Source::addProperty(SourceProperty* property) {
	properties.push_back(property);
}

ref_ptr<Candidate> Source::getCandidate() const {
	ref_ptr<Candidate> candidate = new Candidate();
	for (int i = 0; i < properties.size(); i++)
		(*properties[i]).prepare(*candidate);
	return candidate;
}

void SourceList::addSource(Source* source, double lumi) {
	sources.push_back(source);
	if (luminosities.size() > 0)
		lumi += luminosities.back();
	luminosities.push_back(lumi);
}

ref_ptr<Candidate> SourceList::getCandidate() const {
	if (sources.size() == 0)
		throw std::runtime_error("SourceList: no sources set");
	double r = Random::instance().rand() * luminosities.back();
	int i = 0;
	while ((r > luminosities[i]) and (i < luminosities.size()))
		i++;
	return (sources[i])->getCandidate();
}

SourceParticleType::SourceParticleType(int id) :
		id(id) {
}

void SourceParticleType::prepare(ParticleState& particle) const {
	particle.setId(id);
}

SourceEnergy::SourceEnergy(double energy) :
		E(energy) {
}

void SourceEnergy::prepare(ParticleState& p) const {
	p.setEnergy(E);
}

SourcePowerLawSpectrum::SourcePowerLawSpectrum(double Emin, double Emax,
		double index) :
		Emin(Emin), Emax(Emax), index(index) {
}

void SourcePowerLawSpectrum::prepare(ParticleState& particle) const {
	Random &random = Random::instance();
	double E = random.randPowerLaw(index, Emin, Emax);
	particle.setEnergy(E);
}

void SourceMultipleParticleTypes::add(int id, double a) {
	ids.push_back(id);
	if (abundances.size() > 0)
		a += abundances.back();
	abundances.push_back(a);
}

void SourceMultipleParticleTypes::prepare(ParticleState& particle) const {
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
	int A = massNumberFromNucleusId(id);
	double weightedAbundance = a * pow(A, -index - 1);
	abundance.push_back(weightedAbundance);
	probability.push_back(0);
	normalize();
}

void SourceComposition::add(int A, int Z, double a) {
	add(nucleusId(A, Z), a);
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

void SourcePosition::prepare(ParticleState& particle) const {
	particle.setPosition(position);
}

void SourceMultiplePositions::add(Vector3d pos, double lumi) {
	positions.push_back(pos);
	if (luminosities.size() > 0)
		lumi += luminosities.back();
	luminosities.push_back(lumi);
}

void SourceMultiplePositions::prepare(ParticleState& particle) const {
	if (positions.size() == 0)
		throw std::runtime_error("SourceMultiplePositions: no position set");
	double r = Random().rand() * luminosities.back();
	int i = 0;
	while ((r > luminosities[i]) and (i < luminosities.size()))
		i++;
	particle.setPosition(positions[i]);
}

SourceUniformSphere::SourceUniformSphere(Vector3d center, double radius) :
		center(center), radius(radius) {
}

void SourceUniformSphere::prepare(ParticleState& particle) const {
	Random &random = Random::instance();
	double r = pow(random.rand(), 1. / 3.) * radius;
	particle.setPosition(random.randVector() * r);
}

SourceUniformShell::SourceUniformShell(Vector3d center, double radius) :
		center(center), radius(radius) {
}

void SourceUniformShell::prepare(ParticleState& particle) const {
	Random &random = Random::instance();
	particle.setPosition(random.randVector() * radius);
}

SourceUniformBox::SourceUniformBox(Vector3d origin, Vector3d size) :
		origin(origin), size(size) {
}

void SourceUniformBox::prepare(ParticleState& particle) const {
	Random &random = Random::instance();
	Vector3d pos(random.rand(), random.rand(), random.rand());
	particle.setPosition(pos * size + origin);
}

SourceUniform1D::SourceUniform1D(double minD, double maxD, bool withCosmology) {
	this->withCosmology = withCosmology;
	if (withCosmology) {
		this->minD = comoving2LightTravelDistance(minD);
		this->maxD = comoving2LightTravelDistance(maxD);
	} else {
		this->minD = minD;
		this->maxD = maxD;
	}
}

void SourceUniform1D::prepare(ParticleState& particle) const {
	Random& random = Random::instance();
	double d = random.rand() * (maxD - minD) + minD;
	if (withCosmology)
		d = lightTravel2ComovingDistance(d);
	particle.setPosition(Vector3d(d, 0, 0));
}

SourceDensityGrid::SourceDensityGrid(ref_ptr<ScalarGrid> grid) :
		grid(grid) {
	float sum = 0;
	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				sum += grid->get(ix, iy, iz);
				grid->get(ix, iy, iz) = sum;
			}
		}
	}
	sumDensity = sum;
}

void SourceDensityGrid::prepare(ParticleState& particle) const {
	Random &random = Random::instance();

	// pick random bin; find bin using STL method
	double r = random.rand(sumDensity);
	std::vector<float> &v = grid->getGrid();
	std::vector<float>::iterator it = lower_bound(v.begin(), v.end(), r);
	int i = it - v.begin();
	Vector3d pos = grid->positionFromIndex(i);

	// draw uniform position within bin
	double dx = random.rand() - 0.5;
	double dy = random.rand() - 0.5;
	double dz = random.rand() - 0.5;
	pos += Vector3d(dx, dy, dz) * grid->getSpacing();

	particle.setPosition(pos);
}

SourceDensityGrid1D::SourceDensityGrid1D(ref_ptr<ScalarGrid> grid) :
		grid(grid) {
	if (grid->getNy() != 1)
		throw std::runtime_error("SourceDensityGrid1D: Ny != 1");
	if (grid->getNz() != 1)
		throw std::runtime_error("SourceDensityGrid1D: Nz != 1");

	float sum = 0;
	for (int ix = 0; ix < grid->getNx(); ix++) {
		sum += grid->get(ix, 0, 0);
		grid->get(ix, 0, 0) = sum;
	}
	sumDensity = sum;
}

void SourceDensityGrid1D::prepare(ParticleState& particle) const {
	Random &random = Random::instance();

	// pick random bin; find bin using STL method
	double r = random.rand(sumDensity);
	std::vector<float> &v = grid->getGrid();
	std::vector<float>::iterator it = lower_bound(v.begin(), v.end(), r);
	int i = it - v.begin();
	Vector3d pos = grid->positionFromIndex(i);

	// draw uniform position within bin
	double dx = random.rand() - 0.5;
	pos += Vector3d(dx, 0, 0) * grid->getSpacing();

	particle.setPosition(pos);
}

void SourceIsotropicEmission::prepare(ParticleState& particle) const {
	Random &random = Random::instance();
	particle.setDirection(random.randVector());
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
	Random &random = Random::instance();
	particle.setDirection(random.randConeVector(direction, aperture));
}

SourceRedshift::SourceRedshift(double z) :
		z(z) {
}

void SourceRedshift::prepare(Candidate& candidate) const {
	candidate.setRedshift(z);
}

SourceUniformRedshift::SourceUniformRedshift(double zmin, double zmax) :
		zmin(zmin), zmax(zmax) {
}

void SourceUniformRedshift::prepare(Candidate& candidate) const {
	double z = Random::instance().randUniform(zmin, zmax);
	candidate.setRedshift(z);
}

void SourceRedshift1D::prepare(Candidate& candidate) const {
	double d = candidate.source.getPosition().getMag();
	double z = comovingDistance2Redshift(d);
	candidate.setRedshift(z);
}

} // namespace crpropa

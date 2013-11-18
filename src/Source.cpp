#include "crpropa/Source.h"
#include "crpropa/Random.h"
#include "crpropa/Cosmology.h"

#include "HepPID/ParticleIDMethods.hh"
#include <stdexcept>

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

void SourceList::addSource(Source* source, double weight) {
	sources.push_back(source);
	if (cdf.size() > 0)
		weight += cdf.back();
	cdf.push_back(weight);
}

ref_ptr<Candidate> SourceList::getCandidate() const {
	if (sources.size() == 0)
		throw std::runtime_error("SourceList: no sources set");
	size_t i = Random::instance().randBin(cdf);
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
	particleTypes.push_back(id);
	if (cdf.size() > 0)
		a += cdf.back();
	cdf.push_back(a);
}

void SourceMultipleParticleTypes::prepare(ParticleState& particle) const {
	if (particleTypes.size() == 0)
		throw std::runtime_error("SourceMultipleParticleTypes: no nuclei set");
	size_t i = Random::instance().randBin(cdf);
	particle.setId(particleTypes[i]);
}

SourceComposition::SourceComposition(double Emin, double Rmax, double index) :
		Emin(Emin), Rmax(Rmax), index(index) {
}

void SourceComposition::add(int id, double weight) {
	nuclei.push_back(id);
	int A = massNumber(id);
	int Z = chargeNumber(id);

	double a = 1 + index;
	if (std::abs(a) < std::numeric_limits<double>::min())
		weight *= log(Z * Rmax / Emin);
	else
		weight *= (pow(Z * Rmax, a) - pow(Emin, a)) / a;

	weight *= pow(A, -a);

	if (cdf.size() > 0)
		weight += cdf.back();
	cdf.push_back(weight);
}

void SourceComposition::add(int A, int Z, double a) {
	add(nucleusId(A, Z), a);
}

void SourceComposition::prepare(ParticleState& particle) const {
	if (nuclei.size() == 0)
		throw std::runtime_error("SourceComposition: No source isotope set");

	Random &random = Random::instance();

	// draw random particle type
	size_t i = random.randBin(cdf);
	int id = nuclei[i];
	particle.setId(id);

	// random energy from power law
	int Z = chargeNumber(id);
	particle.setEnergy(random.randPowerLaw(index, Emin, Z * Rmax));
}

SourcePosition::SourcePosition(Vector3d position) :
		position(position) {
}

void SourcePosition::prepare(ParticleState& particle) const {
	particle.setPosition(position);
}

void SourceMultiplePositions::add(Vector3d pos, double weight) {
	positions.push_back(pos);
	if (cdf.size() > 0)
		weight += cdf.back();
	cdf.push_back(weight);
}

void SourceMultiplePositions::prepare(ParticleState& particle) const {
	if (positions.size() == 0)
		throw std::runtime_error("SourceMultiplePositions: no position set");
	size_t i = Random::instance().randBin(cdf);
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
}

void SourceDensityGrid::prepare(ParticleState& particle) const {
	Random &random = Random::instance();

	// draw random bin
	size_t i = random.randBin(grid->getGrid());
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
}

void SourceDensityGrid1D::prepare(ParticleState& particle) const {
	Random &random = Random::instance();

	// draw random bin
	size_t i = random.randBin(grid->getGrid());
	Vector3d pos = grid->positionFromIndex(i);

	// draw uniform position within bin
	double dx = random.rand() - 0.5;
	pos.x += dx * grid->getSpacing();

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
	double d = candidate.source.getPosition().getR();
	double z = comovingDistance2Redshift(d);
	candidate.setRedshift(z);
}

} // namespace crpropa

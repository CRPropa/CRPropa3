#include "crpropa/Source.h"
#include "crpropa/Random.h"
#include "crpropa/Cosmology.h"

#include "HepPID/ParticleIDMethods.hh"

#ifdef CRPROPA_HAVE_MUPARSER
#include "muParser.h"
#endif

#include <stdexcept>

namespace crpropa {

void SourceFeature::prepareParticle(ParticleState& particle) const {
}

void SourceFeature::prepareCandidate(Candidate& candidate) const {
	ParticleState &source = candidate.source;
	prepareParticle(source);
	candidate.created = source;
	candidate.current = source;
	candidate.previous = source;
}

void Source::add(SourceFeature* property) {
	features.push_back(property);
}

ref_ptr<Candidate> Source::getCandidate() const {
	ref_ptr<Candidate> candidate = new Candidate();
	for (int i = 0; i < features.size(); i++)
		(*features[i]).prepareCandidate(*candidate);
	return candidate;
}

void SourceList::add(Source* source, double weight) {
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

void SourceParticleType::prepareParticle(ParticleState& particle) const {
	particle.setId(id);
}

SourceEnergy::SourceEnergy(double energy) :
		E(energy) {
}

void SourceEnergy::prepareParticle(ParticleState& p) const {
	p.setEnergy(E);
}

#ifdef CRPROPA_HAVE_MUPARSER
SourceGenericComposition::SourceGenericComposition(double Emin, double Emax, std::string expression, size_t steps) :
	Emin(Emin), Emax(Emax), expression(expression), steps(steps) {

	// precalculate energy bins
	double logEmin = ::log10(Emin);
	double logEmax = ::log10(Emax);
	double logStep = (logEmax - logEmin) / (steps-1);
	energy.resize(steps);
	for (size_t i = 0; i < steps; i++) {
		energy[i] = ::pow(10, logEmin + i * logStep);
	}
}

void SourceGenericComposition::add(int id, double weight) {
	int A = massNumber(id);
	int Z = chargeNumber(id);

	Nucleus n;
	n.id = id;

	// calculate nuclei cdf
	mu::Parser p;
	double E;
    p.DefineVar("E", &E); 
    p.DefineConst("Emin", Emin); 
    p.DefineConst("Emax", Emax); 
    p.DefineConst("steps", steps); 
    p.DefineConst("A", (double)A); 
    p.DefineConst("Z", (double)Z); 
    p.SetExpr(expression);

	// calculate pdf
	n.cdf.resize(steps);
    for (std::size_t i=0; i<steps; ++i) {
		E = energy[i];
		n.cdf[i] = p.Eval();
    }

	// integrate
    for (std::size_t i=1; i<steps; ++i) {
		n.cdf[i] = (n.cdf[i-1] + n.cdf[i]) * (energy[i] - energy[i-1]) / 2;
    }

	// cumulate
	n.cdf[0] = 0;
    for (std::size_t i=1; i<steps; ++i) {
		n.cdf[i] += n.cdf[i-1];
    }

	nuclei.push_back(n);

	if (cdf.size() > 0)
		weight += cdf.back();
	cdf.push_back(weight);
}

void SourceGenericComposition::add(int A, int Z, double a) {
	add(nucleusId(A, Z), a);
}

void SourceGenericComposition::prepareParticle(ParticleState& particle) const {
	if (nuclei.size() == 0)
		throw std::runtime_error("SourceComposition: No source isotope set");

	Random &random = Random::instance();


	// draw random particle type
	size_t iN = random.randBin(cdf);
	const Nucleus &n = nuclei[iN];
	particle.setId(n.id);

	// random energy
	double E = interpolate(random.rand() * n.cdf.back(), n.cdf, energy);
	particle.setEnergy(E);
}
#endif

SourcePowerLawSpectrum::SourcePowerLawSpectrum(double Emin, double Emax,
		double index) :
		Emin(Emin), Emax(Emax), index(index) {
}

void SourcePowerLawSpectrum::prepareParticle(ParticleState& particle) const {
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

void SourceMultipleParticleTypes::prepareParticle(ParticleState& particle) const {
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

void SourceComposition::prepareParticle(ParticleState& particle) const {
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

void SourcePosition::prepareParticle(ParticleState& particle) const {
	particle.setPosition(position);
}

void SourceMultiplePositions::add(Vector3d pos, double weight) {
	positions.push_back(pos);
	if (cdf.size() > 0)
		weight += cdf.back();
	cdf.push_back(weight);
}

void SourceMultiplePositions::prepareParticle(ParticleState& particle) const {
	if (positions.size() == 0)
		throw std::runtime_error("SourceMultiplePositions: no position set");
	size_t i = Random::instance().randBin(cdf);
	particle.setPosition(positions[i]);
}

SourceUniformSphere::SourceUniformSphere(Vector3d center, double radius) :
		center(center), radius(radius) {
}

void SourceUniformSphere::prepareParticle(ParticleState& particle) const {
	Random &random = Random::instance();
	double r = pow(random.rand(), 1. / 3.) * radius;
	particle.setPosition(random.randVector() * r);
}

SourceUniformShell::SourceUniformShell(Vector3d center, double radius) :
		center(center), radius(radius) {
}

void SourceUniformShell::prepareParticle(ParticleState& particle) const {
	Random &random = Random::instance();
	particle.setPosition(random.randVector() * radius);
}

SourceUniformBox::SourceUniformBox(Vector3d origin, Vector3d size) :
		origin(origin), size(size) {
}

void SourceUniformBox::prepareParticle(ParticleState& particle) const {
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

void SourceUniform1D::prepareParticle(ParticleState& particle) const {
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

void SourceDensityGrid::prepareParticle(ParticleState& particle) const {
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

void SourceDensityGrid1D::prepareParticle(ParticleState& particle) const {
	Random &random = Random::instance();

	// draw random bin
	size_t i = random.randBin(grid->getGrid());
	Vector3d pos = grid->positionFromIndex(i);

	// draw uniform position within bin
	double dx = random.rand() - 0.5;
	pos.x += dx * grid->getSpacing();

	particle.setPosition(pos);
}

void SourceIsotropicEmission::prepareParticle(ParticleState& particle) const {
	Random &random = Random::instance();
	particle.setDirection(random.randVector());
}

SourceDirection::SourceDirection(Vector3d direction) :
		direction(direction) {
}

void SourceDirection::prepareParticle(ParticleState& particle) const {
	particle.setDirection(direction);
}

SourceEmissionCone::SourceEmissionCone(Vector3d direction, double aperture) :
		direction(direction), aperture(aperture) {
}

void SourceEmissionCone::prepareParticle(ParticleState& particle) const {
	Random &random = Random::instance();
	particle.setDirection(random.randConeVector(direction, aperture));
}

SourceRedshift::SourceRedshift(double z) :
		z(z) {
}

void SourceRedshift::prepareCandidate(Candidate& candidate) const {
	candidate.setRedshift(z);
}

SourceUniformRedshift::SourceUniformRedshift(double zmin, double zmax) :
		zmin(zmin), zmax(zmax) {
}

void SourceUniformRedshift::prepareCandidate(Candidate& candidate) const {
	double z = Random::instance().randUniform(zmin, zmax);
	candidate.setRedshift(z);
}

void SourceRedshift1D::prepareCandidate(Candidate& candidate) const {
	double d = candidate.source.getPosition().getR();
	double z = comovingDistance2Redshift(d);
	candidate.setRedshift(z);
}

} // namespace crpropa

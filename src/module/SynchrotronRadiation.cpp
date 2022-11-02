#include "crpropa/module/SynchrotronRadiation.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

SynchrotronRadiation::SynchrotronRadiation(ref_ptr<MagneticField> field, bool havePhotons, double thinning, int nSamples, double limit) {
	setField(field);
	setBrms(0);
	initSpectrum();
	setHavePhotons(havePhotons);
	setLimit(limit);
	setSecondaryThreshold(1e6 * eV);
	setMaximumSamples(nSamples);
}

SynchrotronRadiation::SynchrotronRadiation(double Brms, bool havePhotons, double thinning, int nSamples, double limit) {
	setBrms(Brms);
	initSpectrum();
	setHavePhotons(havePhotons);
	setLimit(limit);
	setSecondaryThreshold(1e6 * eV);
	setMaximumSamples(nSamples);
}

void SynchrotronRadiation::setField(ref_ptr<MagneticField> f) {
	this->field = f;
}

ref_ptr<MagneticField> SynchrotronRadiation::getField() {
	return field;
}

void SynchrotronRadiation::setBrms(double Brms) {
	this->Brms = Brms;
}

double SynchrotronRadiation::getBrms() {
	return Brms;
}

void SynchrotronRadiation::setHavePhotons(bool havePhotons) {
	this->havePhotons = havePhotons;
}

bool SynchrotronRadiation::getHavePhotons() {
	return havePhotons;
}

void SynchrotronRadiation::setThinning(double thinning) {
	this->thinning = thinning;
}

double SynchrotronRadiation::getThinning() {
	return thinning;
}

void SynchrotronRadiation::setLimit(double limit) {
	this->limit = limit;
}

double SynchrotronRadiation::getLimit() {
	return limit;
}

void SynchrotronRadiation::setMaximumSamples(int nmax) {
	maximumSamples = nmax;
}

int SynchrotronRadiation::getMaximumSamples() {
	return maximumSamples;
}

void SynchrotronRadiation::setSecondaryThreshold(double threshold) {
	secondaryThreshold = threshold;
}

double SynchrotronRadiation::getSecondaryThreshold() const {
	return secondaryThreshold;
}

void SynchrotronRadiation::initSpectrum() {
	std::string filename = getDataPath("Synchrotron/spectrum.txt");
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error("SynchrotronRadiation: could not open file " + filename);

	// clear previously loaded interaction rates
	tabx.clear();
	tabCDF.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabx.push_back(pow(10, a));
				tabCDF.push_back(b);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
}

void SynchrotronRadiation::process(Candidate *candidate) const {
	double charge = fabs(candidate->current.getCharge());
	if (charge == 0)
		return; // only charged particles

	// calculate gyroradius, evaluated at the current position
	double z = candidate->getRedshift();
	double B;
	if (field.valid()) {
		Vector3d Bvec = field->getField(candidate->current.getPosition(), z);
		B = Bvec.cross(candidate->current.getDirection()).getR();
	} else {
		B = sqrt(2. / 3) * Brms; // average perpendicular field component
	}
	B *= pow(1 + z, 2); // cosmological scaling
	double Rg = candidate->current.getMomentum().getR() / charge / B;

	// calculate energy loss
	double lf = candidate->current.getLorentzFactor();
	double dEdx = 1. / 6 / M_PI / epsilon0 * pow(lf * lf - 1, 2) * pow(charge / Rg, 2); // Jackson p. 770 (14.31)
	double step = candidate->getCurrentStep() / (1 + z); // step size in local frame
	double dE = step * dEdx;

	// apply energy loss and limit next step
	double E = candidate->current.getEnergy();
	candidate->current.setEnergy(E - dE);
	candidate->limitNextStep(limit * E / dEdx);

	// optionally add secondary photons
	if (not(havePhotons))
		return;

	// check if photons with energies > 14 * Ecrit are possible
	double Ecrit = 3. / 4 * h_planck / M_PI * c_light * pow(lf, 3) / Rg;
	if (14 * Ecrit < secondaryThreshold)
		return;

	// draw photons up to the total energy loss
	// if maximumSamples is reached before that, compensate the total energy afterwards
	Random &random = Random::instance();
	double dE0 = dE;
	std::vector<double> energies;
	int counter = 0;
	while (dE > 0) {
		// draw random value between 0 and maximum of corresponding cdf
		// choose bin of s where cdf(x) = cdf_rand -> x_rand
		size_t i = random.randBin(tabCDF); // draw random bin (upper bin boundary returned)
		double binWidth = (tabx[i] - tabx[i-1]);
		double x = tabx[i-1] + random.rand() * binWidth; // draw random x uniformly distributed in bin
		double Ephoton = x * Ecrit;

		// if the remaining energy is not sufficient check for random accepting
		if (Ephoton > dE) {
			if (random.rand() > (dE / Ephoton))
				break; // not accepted
		}

		// only activate the "per-step" sampling if maximumSamples is explicitly set.
		if (maximumSamples > 0) {
			if (counter >= maximumSamples) 
				break;			
		}

		// store energies in array
		energies.push_back(Ephoton);

		// energy loss
		dE -= Ephoton;

		// counter for sampling break condition;
		counter++;
	}

	// while loop before gave total energy which is just a fraction of the required
	double w1 = 1;
	if (maximumSamples > 0 && dE > 0)
		w1 = 1. / (1. - dE / dE0); 

	// loop over sampled photons and attribute weights accordingly
	for (int i = 0; i < energies.size(); i++) {
		double Ephoton = energies[i];
		double f = Ephoton / (E - dE0);
		double w = w1 / pow(f, thinning);

		// thinning procedure: accepts only a few random secondaries
		if (random.rand() < pow(f, thinning)) {
			Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
			if (Ephoton > secondaryThreshold) // create only photons with energies above threshold
				candidate->addSecondary(22, Ephoton, pos, w, interactionTag);
		}
	}
}

std::string SynchrotronRadiation::getDescription() const {
	std::stringstream s;
	s << "Synchrotron radiation";
	if (field.valid())
		s << " for specified magnetic field";
	else
		s << " for Brms = " << Brms / nG << " nG";
	if (havePhotons)
		s << ", synchrotron photons E > " << secondaryThreshold / eV << " eV";
	else
		s << ", no synchrotron photons";
	if (maximumSamples > 0)
		s << "maximum number of photon samples: " << maximumSamples;
	if (thinning > 0)
		s << "thinning parameter: " << thinning; 
	return s.str();
}

void SynchrotronRadiation::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string SynchrotronRadiation::getInteractionTag() const {
	return interactionTag;
}

} // namespace crpropa

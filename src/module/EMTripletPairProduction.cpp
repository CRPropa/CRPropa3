#include "crpropa/module/EMTripletPairProduction.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

static const double mec2 = mass_electron * c_squared;

EMTripletPairProduction::EMTripletPairProduction(PhotonField photonField, bool haveElectrons, double thinning, double limit) {
	setPhotonField(photonField);
	this->haveElectrons = haveElectrons;
	this->limit = limit;
	this->thinning = thinning;
}

void EMTripletPairProduction::setPhotonField(PhotonField photonField) {
	this->photonField = photonField;
	std::string fname = photonFieldName(photonField);
	setDescription("EMTripletPairProduction: " + fname);
	initRate(getDataPath("EMTripletPairProduction/rate_" + fname + ".txt"));
	initCumulativeRate(getDataPath("EMTripletPairProduction/cdf_" + fname + ".txt"));
}

void EMTripletPairProduction::setHaveElectrons(bool haveElectrons) {
	this->haveElectrons = haveElectrons;
}

void EMTripletPairProduction::setLimit(double limit) {
	this->limit = limit;
}

void EMTripletPairProduction::setThinning(double thinning) {
	this->thinning = thinning;
}

void EMTripletPairProduction::initRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error("EMTripletPairProduction: could not open file " + filename);

	// clear previously loaded interaction rates
	tabEnergy.clear();
	tabRate.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabEnergy.push_back(pow(10, a) * eV);
				tabRate.push_back(b / Mpc);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
}

void EMTripletPairProduction::initCumulativeRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error(
				"EMTripletPairProduction: could not open file " + filename);

	// clear previously loaded tables
	tabE.clear();
	tabs.clear();
	tabCDF.clear();
	
	// skip header
	while (infile.peek() == '#')
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');

	// read s values in first line
	double a;
	infile >> a; // skip first value
	while (infile.good() and (infile.peek() != '\n')) {
		infile >> a;
		tabs.push_back(pow(10, a) * eV * eV);
	}

	// read all following lines: E, cdf values
	while (infile.good()) {
		infile >> a;
		if (!infile)
			break;  // end of file
		tabE.push_back(pow(10, a) * eV);
		std::vector<double> cdf;
		for (int i = 0; i < tabs.size(); i++) {
			infile >> a;
			cdf.push_back(a / Mpc);
		}
		tabCDF.push_back(cdf);
	}
	infile.close();
}

void EMTripletPairProduction::performInteraction(Candidate *candidate) const {

	if  (abs(candidate->current.getId()) != 11)
		return;

	// scale the particle energy instead of background photons
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	if (E < tabE.front() or E > tabE.back())
		return;

	// sample the value of eps
	Random &random = Random::instance();
	size_t i = closestIndex(E, tabE);
	size_t j = random.randBin(tabCDF[i]);
	double s_kin = pow(10, log10(tabs[j]) + (random.rand() - 0.5) * 0.1);
	double eps = s_kin / 4 / E; // random background photon energy

	// Use approximation from A. Mastichiadis et al., Astroph. Journ. 300:178-189 (1986), eq. 30.
	// This approx is valid only for alpha >=100 where alpha = p0*eps*costheta - E0*eps
	// For our purposes, me << E0 --> p0~E0 --> alpha = E0*eps*(costheta - 1) >= 100
	double Epp = 5.7e-1 * pow(eps / mec2, -0.56) * pow(E / mec2, 0.44) * mec2;

	double w0 = candidate->getWeight();
	double f = Epp / E;

	if (haveElectrons) {
		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
		if (random.rand() < pow(1 - f, thinning)) {
			double w = w0 / pow(1 - f, thinning);
			candidate->addSecondary(11, Epp, pos, w);
		}
		if (random.rand() < pow(f, thinning)) {
			double w = w0 / pow(f, thinning);
			candidate->addSecondary(-11, Epp, pos, w);
		}
	}

	// update the primary particle energy; do this after adding the secondaries to correctly set the secondaries parent
	candidate->current.setEnergy((E - 2 * Epp));
}

void EMTripletPairProduction::process(Candidate *candidate) const {
	// check if electron / positron
	int id = candidate->current.getId();
	if (abs(id) != 11)
		return;

	// scale the particle energy instead of background photons
	double z = candidate->getRedshift();
	double E = (1 + z) * candidate->current.getEnergy();

	// check if in tabulated energy range
	if (E < tabEnergy.front() or (E > tabEnergy.back()))
		return;

	// cosmological scaling of interaction distance (comoving)
	double scaling = pow(1 + z, 2) * photonFieldScaling(photonField, z);
	double rate = scaling * interpolate(E, tabEnergy, tabRate);

	// run this loop at least once to limit the step size
	double step = candidate->getCurrentStep();
	while (step > 0) {
		// check for interaction
		Random &random = Random::instance();
		double randDistance = -log(random.rand()) / rate;
		if (step < randDistance) {
			candidate->limitNextStep(limit / rate);
			return;
		}
		performInteraction(candidate);

		step -= randDistance;
	}
}

} // namespace crpropa

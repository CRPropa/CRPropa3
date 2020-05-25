#include "crpropa/module/EMDoublePairProduction.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

EMDoublePairProduction::EMDoublePairProduction(PhotonField photonField, bool haveElectrons, double thinning, double limit) {
	setPhotonField(photonField);
	this->haveElectrons = haveElectrons;
	this->thinning = thinning;
	this->limit = limit;
}

void EMDoublePairProduction::setPhotonField(PhotonField photonField) {
	this->photonField = photonField;
	std::string fname = photonFieldName(photonField);
	setDescription("EMDoublePairProduction: " + fname);
	initRate(getDataPath("EMDoublePairProduction/rate_" + fname + ".txt"));
}

void EMDoublePairProduction::setHaveElectrons(bool haveElectrons) {
	this->haveElectrons = haveElectrons;
}

void EMDoublePairProduction::setLimit(double limit) {
	this->limit = limit;
}

void EMDoublePairProduction::setThinning(double thinning) {
	this->thinning = thinning;
}

void EMDoublePairProduction::initRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error("EMDoublePairProduction: could not open file " + filename);

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


void EMDoublePairProduction::performInteraction(Candidate *candidate) const {
	// the photon is lost in interaction
	candidate->setActive(false);

	if (not haveElectrons)
		return;

	// Use assumption of Lee 96 arXiv:9604098
	// Energy is equally shared between one e+e- pair, but take mass of second e+e- pair into account.
	// This approximation has been shown to be valid within -1.5%.
	double E = candidate->current.getEnergy();
	double Ee = (E - 2 * mass_electron * c_squared) / 2;

	Random &random = Random::instance();
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());

	double f = Ee / E;
	double w0 = candidate->getWeight();

	if (haveElectrons) {
		if (random.rand() < pow(1 - f, thinning)) {
			double w = w0 / pow(1 - f, thinning);
			candidate->addSecondary( 11, Ee, pos, w);
		} 
		if (random.rand() < pow(f, thinning)) {
			double w = w0 / pow(f, thinning);
			candidate->addSecondary(-11, Ee, pos, w);
		}
	}
}

void EMDoublePairProduction::process(Candidate *candidate) const {
	// check if photon
	if (candidate->current.getId() != 22)
		return;

	// scale the electron energy instead of background photons
	double z = candidate->getRedshift();
	double E = (1 + z) * candidate->current.getEnergy();

	// check if in tabulated energy range
	if (E < tabEnergy.front() or (E > tabEnergy.back()))
		return;

	// interaction rate
	double rate = interpolate(E, tabEnergy, tabRate);
	rate *= pow(1 + z, 2) * photonFieldScaling(photonField, z);

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

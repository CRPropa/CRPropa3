#include "crpropa/module/EMDoublePairProduction.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

EMDoublePairProduction::EMDoublePairProduction(PhotonField photonField,
											   bool haveElectrons,
											   std::string tag,
											   double limit) {
	setPhotonField(photonField);
	this->spaceTimeGrid = ScalarGrid4d();
	this->spaceGrid = ScalarGrid();
	this->haveElectrons = haveElectrons;
	this->tag = tag;
	this->limit = limit;
	setDescription("EMDoublePairProduction_isotropicConstant");
}

EMDoublePairProduction::EMDoublePairProduction(PhotonField photonField,
											   ScalarGrid4d spaceTimeGrid,
											   bool haveElectrons,
											   std::string tag,
											   double limit) {
	setPhotonField(photonField);
	this->spaceTimeGrid = spaceTimeGrid;
	this->spaceGrid = ScalarGrid();
	this->haveElectrons = haveElectrons;
	this->tag = tag;
	this->limit = limit;
	setDescription("EMDoublePairProduction_spaceTimeDependent");
}

EMDoublePairProduction::EMDoublePairProduction(PhotonField photonField,
											   ScalarGrid spaceGrid,
											   bool haveElectrons,
											   std::string tag,
											   double limit) {
	setPhotonField(photonField);
	this->spaceTimeGrid = ScalarGrid4d();
	this->spaceGrid = spaceGrid;
	this->haveElectrons = haveElectrons;
	this->tag = tag;
	this->limit = limit;
	setDescription("EMDoublePairProduction_spaceDependentConstant");
	
}

void EMDoublePairProduction::setPhotonField(PhotonField photonField) {
	this->photonField = photonField;
	std::string fname = photonFieldName(photonField);
	initRate(getDataPath("EMDoublePairProduction/rate_" + fname + ".txt"));
}

void EMDoublePairProduction::setHaveElectrons(bool haveElectrons) {
	this->haveElectrons = haveElectrons;
}

void EMDoublePairProduction::setLimit(double limit) {
	this->limit = limit;
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

	candidate->addSecondary( 11, Ee, pos, tag);
	candidate->addSecondary(-11, Ee, pos, tag);
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
	double rate = 1.;

	// geometric scaling
	Vector3d pos = candidate->current.getPosition();
	double time = candidate->getTrajectoryLength()/c_light;
	
	const std::string description = getDescription();
	if (description == "EMDoublePairProduction_isotropicConstant") {
		// do nothing, just check for correct initialization
	} else if (description == "EMDoublePairProduction_spaceDependentConstant") {
		rate *= spaceGrid.interpolate(pos);
	} else if (description == "EMDoublePairProduction_spaceTimeDependent") {
		rate *= spaceTimeGrid.interpolate(pos, time);
	} else {
		throw std::runtime_error("EMDoublePairProduction: invalid description string");
	}
	if (rate == 0.)
		return;

	rate *= interpolate(E, tabEnergy, tabRate);
	rate *= pow(1 + z, 2) * photonFieldScaling(photonField, z);

	// check for interaction
	Random &random = Random::instance();
	double randDistance = -log(random.rand()) / rate;
	if (candidate->getCurrentStep() > randDistance)
		performInteraction(candidate);
	else
		candidate->limitNextStep(limit / rate);
}

} // namespace crpropa

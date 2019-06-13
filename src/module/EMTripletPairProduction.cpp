#include "crpropa/module/EMTripletPairProduction.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

static const double mec2 = mass_electron * c_squared;

EMTripletPairProduction::EMTripletPairProduction(PhotonField photonField,
												 bool haveElectrons,
												 std::string tag,
												 double limit) {
	setPhotonField(photonField);
	this->spaceTimeGrid = ScalarGrid4d();
	this->spaceGrid = ScalarGrid();
	this->haveElectrons = haveElectrons;
	this->tag = tag;
	this->limit = limit;
	setDescription("EMTripletPairProduction_isotropicConstant");
}

EMTripletPairProduction::EMTripletPairProduction(PhotonField photonField,
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
	setDescription("EMTripletPairProduction_spaceTimeDependent");
}

EMTripletPairProduction::EMTripletPairProduction(PhotonField photonField,
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
	setDescription("EMTripletPairProduction_spaceDependentConstant");
}

void EMTripletPairProduction::setPhotonField(PhotonField photonField) {
	this->photonField = photonField;
	std::string fname = photonFieldName(photonField);
	initRate(getDataPath("EMTripletPairProduction/rate_" + fname + ".txt"));
	initCumulativeRate(getDataPath("EMTripletPairProduction/cdf_" + fname + ".txt"));
}

void EMTripletPairProduction::setHaveElectrons(bool haveElectrons) {
	this->haveElectrons = haveElectrons;
}

void EMTripletPairProduction::setLimit(double limit) {
	this->limit = limit;
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
	// scale the particle energy instead of background photons
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	if (E < tabE.front() or E > tabE.back())
		return;

	// geometric scaling
		// geometric scaling
	Vector3d pos = candidate->current.getPosition();
	const double time = candidate->getTrajectoryLength()/c_light;
	
	double geometricScaling = 1.;
	const std::string description = getDescription();
	if (description == "EMPairProduction_isotropicConstant") {
		// do nothing, just check for correct initialization
	} else if (description == "EMPairProduction_spaceDependentConstant") {
		geometricScaling *= spaceGrid.interpolate(pos);
	} else if (description == "EMPairProduction_spaceTimeDependent") {
		geometricScaling *= spaceTimeGrid.interpolate(pos, time);
	} else {
		throw std::runtime_error("EMPairProduction: invalid description string");
	}
	if (geometricScaling == 0.)
		return;

	// sample the value of eps
	Random &random = Random::instance();
	size_t i = closestIndex(E, tabE);
	size_t j = random.randBin(tabCDF[i]);
	double s_kin = pow(10, log10(tabs[j] * geometricScaling) + (random.rand() - 0.5) * 0.1);
	double eps = s_kin / 4 / E; // random background photon energy

	// Use approximation from A. Mastichiadis et al., Astroph. Journ. 300:178-189 (1986), eq. 30.
	// This approx is valid only for alpha >=100 where alpha = p0*eps*costheta - E0*eps
	// For our purposes, me << E0 --> p0~E0 --> alpha = E0*eps*(costheta - 1) >= 100  <= Note by Mario: How is this even supposed to be >0?
	double Epp = 5.7e-1 * pow(eps/mec2, -0.56) * pow(E/mec2, 0.44) * mec2;

	if (haveElectrons) {
		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
		candidate->addSecondary( 11, Epp, pos, tag);
		candidate->addSecondary(-11, Epp, pos, tag);
	}

	// update the primary particle energy; do this after adding the secondaries to correctly set the secondaries parent
	candidate->current.setEnergy(E - 2 * Epp);
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

	// geometric scaling
	Vector3d pos = candidate->current.getPosition();
    double time = candidate->getTrajectoryLength()/c_light;
	double rate = spaceTimeGrid.interpolate(pos, time);
	if (rate == 0.)
		return;
	
	rate *= interpolate(E, tabEnergy, tabRate);
	// cosmological scaling of interaction distance (comoving)
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

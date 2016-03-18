#include "crpropa/module/EMDoublePairProduction.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

EMDoublePairProduction::EMDoublePairProduction(PhotonField photonField,
		bool haveElectrons, double limit) {
	setPhotonField(photonField);
	this->haveElectrons = haveElectrons;
	this->limit = limit;
}

void EMDoublePairProduction::setPhotonField(PhotonField photonField) {
	this->photonField = photonField;
	switch (photonField) {
	case CMB:
		setDescription("EMDoublePairProduction: CMB");
		initRate(getDataPath("EMDoublePairProduction_CMB.txt"));
		break;
	case IRB:  // default: Kneiske '04 IRB model
	case IRB_Kneiske04:
		setDescription("EMDoublePairProduction: IRB (Kneiske 2004)");
		initRate(getDataPath("EMDoublePairProduction_IRB_Kneiske04.txt"));
		break;
	case IRB_Stecker05:
		setDescription("EMDoublePairProduction: IRB (Stecker 2005)");
		initRate(getDataPath("EMDoublePairProduction_IRB_Stecker05.txt"));
		break;
	case IRB_Franceschini08:
		setDescription("EMDoublePairProduction: IRB (Franceschini 2008)");
		initRate(getDataPath("EMDoublePairProduction_IRB_Franceschini08.txt"));
		break;
	case IRB_Finke10:
		setDescription("EMDoublePairProduction: IRB (Finke 2010)");
		initRate(getDataPath("EMDoublePairProduction_IRB_Finke10.txt"));
		break;
	case IRB_Dominguez11:
		setDescription("EMDoublePairProduction: IRB (Dominguez 2011)");
		initRate(getDataPath("EMDoublePairProduction_IRB_Dominguez11.txt"));
		break;
	case IRB_Gilmore12:
		setDescription("EMDoublePairProduction: IRB (Gilmore 2012)");
		initRate(getDataPath("EMDoublePairProduction_IRB_Gilmore12.txt"));
		break;
	case URB_Protheroe96:
		setDescription("EMDoublePairProduction: URB (Protheroe 1996)");
		initRate(getDataPath("EMDoublePairProduction_URB_Protheroe96.txt"));
		break;
	default:
		throw std::runtime_error(
				"EMDoublePairProduction: unknown photon background");
	}
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
		throw std::runtime_error(
				"EMDoublePairProduction: could not open file " + filename);

	// clear previously loaded interaction rates
	tabPhotonEnergy.clear();
	tabInteractionRate.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabPhotonEnergy.push_back(pow(10, a) * eV);
				tabInteractionRate.push_back(b / Mpc);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
}

void EMDoublePairProduction::performInteraction(Candidate *candidate) const {
	double E = candidate->current.getEnergy();

	if (haveElectrons){
		double Ee = (E-2.*mass_electron*c_squared)/2.; // Use assumption of Lee 96 (i.e., all the energy goes equaly shared between only 1 couple of e+e- but take mass of second e+e- pair into account. In DPPpaper has been shown that this approximation is valid within -1.5%
		Vector3d pos = randomPositionInPropagationStep(candidate);
		candidate->addSecondary(11, Ee, pos);
		candidate->addSecondary(-11, Ee, pos);
	}
	candidate->setActive(false);
}

void EMDoublePairProduction::process(Candidate *candidate) const {
	double step = candidate->getCurrentStep();
	double z = candidate->getRedshift();

	// check if photon
	int id = candidate->current.getId();
	if (id != 22)
		return;

	// instead of scaling the background photon energies, scale the photon energy
	double E = (1 + z) * candidate->current.getEnergy();

	// check if in tabulated energy range
	if (E < tabPhotonEnergy.front() or (E > tabPhotonEnergy.back()))
		return;

	// find interaction with minimum random distance
	Random &random = Random::instance();
	double randDistance = std::numeric_limits<double>::max();

	// comological scaling of interaction distance (comoving)
	double scaling = pow(1 + z, 3) * photonFieldScaling(photonField, z);
	double rate = scaling * interpolate(E, tabPhotonEnergy, tabInteractionRate);
	randDistance = -log(random.rand()) / rate;

	candidate->limitNextStep(limit / rate);
	// check if interaction does not happen
	if (step < randDistance) {
		return;
	}

	// interact
	performInteraction(candidate);
}

} // namespace crpropa

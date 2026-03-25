#include "crpropa/module/EMDoublePairProduction.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"
#include "crpropa/Common.h"

#include <fstream>
#include <locale>
#include <limits>
#include <stdexcept>
#include <filesystem>
#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>

#if defined(__APPLE__) && defined(_LIBCPP_VERSION)
	namespace fs = std::__fs::filesystem;
#else
	namespace fs = std::filesystem;
#endif

namespace crpropa {

EMDoublePairProduction::EMDoublePairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons, double thinning, double limit, ref_ptr<Surface> surface) {

	setSurface(surface);
	setPhotonField(photonField);
	setHaveElectrons(haveElectrons);
	setLimit(limit);
	setThinning(thinning);
	
}

void EMDoublePairProduction::setPhotonField(ref_ptr<PhotonField> photonField) {
	
	this->photonField = photonField;
	std::string fname = photonField->getFieldName();
	setDescription("EMDoublePairProduction: " + fname);
	
	// choose the right interaction rates for the used photon field
	if (!this->photonField->hasPositionDependence()) {
		this->interactionRates = new InteractionRatesHomogeneous(
			getDataPath("EMDoublePairProduction/rate_" + fname + ".txt")
		);
		
	} else {
		this->interactionRates = new InteractionRatesPositionDependent(
			getDataPath("EMDoublePairProduction/"+fname+"/Rate/"),
			"",
			this->surface
		);
		
	}
	
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

void EMDoublePairProduction::setSurface(ref_ptr<Surface> surface) {
		this->surface = surface;
}

ref_ptr<Surface> EMDoublePairProduction::getSurface() const {
		return this->surface;
}

void EMDoublePairProduction::setInteractionRates(ref_ptr<InteractionRates> intRates) {
	this->interactionRates = intRates;
}

ref_ptr<InteractionRates> EMDoublePairProduction::getInteractionRates() const {
	return this->interactionRates;
}

void EMDoublePairProduction::initRate(std::string path) {
	this->interactionRates->initRate(path);
}

void EMDoublePairProduction::performInteraction(Candidate *candidate) const {

	// Use assumption of Lee 96 arXiv:9604098
	// Energy is equally shared between one e+e- pair, but take mass of second e+e- pair into account.
	// This approximation has been shown to be valid within -1.5%.
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);
	double Ee = (E - 2 * mass_electron * c_squared) / 2;

	// the photon is lost after the interaction
	candidate->setActive(false);

	if (not haveElectrons)
		return;
	
	Random &random = Random::instance();
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());

	double f = Ee / E;

		if (random.rand() < pow(1 - f, thinning)) {
			double w = 1. / pow(1 - f, thinning);
			candidate->addSecondary( 11, Ee / (1 + z), pos, w, interactionTag);
		} 
		if (random.rand() < pow(f, thinning)) {
			double w = 1. / pow(f, thinning);
			candidate->addSecondary(-11, Ee / (1 + z), pos, w, interactionTag);
		}

}

void EMDoublePairProduction::process(Candidate *candidate) const {
	
	// check if photon
	if (candidate->current.getId() != 22)
		return;

	// scale the electron energy instead of background photons
	double z = candidate->getRedshift();
	double E = (1 + z) * candidate->current.getEnergy();
	Vector3d position = candidate->current.getPosition();

	// interaction rate
	double rate = this->interactionRates->getProcessRate(E, position);
		
	if (rate < 0)
		return;
		
	rate *= pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);

	// check for interaction
	Random &random = Random::instance();
	double randDistance = -log(random.rand()) / rate;
	double step = candidate->getCurrentStep();
	if (step < randDistance) {
		candidate->limitNextStep(limit / rate);
		return;
	} else { // after performing interaction photon ceases to exist (hence return)
		performInteraction(candidate);
		return;
	}

}

void EMDoublePairProduction::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string EMDoublePairProduction::getInteractionTag() const {
	return interactionTag;
}


} // namespace crpropa

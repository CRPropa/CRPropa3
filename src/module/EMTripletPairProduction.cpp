#include "crpropa/module/EMTripletPairProduction.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"
#include "crpropa/Common.h"

#include <fstream>
#include <locale>
#include <limits>
#include <stdexcept>
#include <filesystem>

#if defined(__APPLE__) && defined(_LIBCPP_VERSION)
	namespace fs = std::__fs::filesystem;
#else
	namespace fs = std::filesystem;
#endif

namespace crpropa {

static const double mec2 = mass_electron * c_squared;

EMTripletPairProduction::EMTripletPairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons, double thinning, double limit, ref_ptr<Surface> surface) {

	setSurface(surface);
	setPhotonField(photonField);
	setHaveElectrons(haveElectrons);
	setLimit(limit);
	setThinning(thinning);
	
}

void EMTripletPairProduction::setPhotonField(ref_ptr<PhotonField> photonField) {
	this->photonField = photonField;
	std::string fname = photonField->getFieldName();
	setDescription("EMTripletPairProduction: " + fname);
	if (!this->photonField->hasPositionDependence()){
		this->interactionRates = new InteractionRatesHomogeneous(
			getDataPath("EMTripletPairProduction/rate_" + fname + ".txt"),
			getDataPath("EMTripletPairProduction/cdf_" + fname + ".txt")
		);

	} else {
		this->interactionRates = new InteractionRatesPositionDependent(
			getDataPath("EMTripletPairProduction/"+fname+"/Rate/"),
			getDataPath("EMTripletPairProduction/"+fname+"/CumulativeRate/"),
			this->surface
		);		
	}
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

void EMTripletPairProduction::setSurface(ref_ptr<Surface> surface) {
		this->surface = surface;
}

ref_ptr<Surface> EMTripletPairProduction::getSurface() const {
		return this->surface;
}

void EMTripletPairProduction::setInteractionRates(ref_ptr<InteractionRates> intRates) {
	this->interactionRates = intRates;
}

ref_ptr<InteractionRates> EMTripletPairProduction::getInteractionRates() const {
	return this->interactionRates;
}

void EMTripletPairProduction::initRate(std::string path) {
	this->interactionRates->initRate(path);
}

void EMTripletPairProduction::initCumulativeRate(std::string path) {
	this->interactionRates->initCumulativeRate(path);
}

void EMTripletPairProduction::performInteraction(Candidate *candidate) const {
	
	// scale the particle energy instead of background photons
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);
	Vector3d position = candidate->current.getPosition();
	
	std::vector<double> tabE;
	std::vector<double> tabs;
	std::vector<std::vector<double>> tabCDF;
	
	this->interactionRates->loadPerformInteractionTabs(position, tabE, tabs, tabCDF);
	
	if (E < tabE.front() or E > tabE.back())
		return;
	
	// sample the value of eps
	Random &random = Random::instance();
	size_t i = closestIndex(E, tabE);
	size_t j = random.randBin(tabCDF[i]);
	double s_kin = pow(10, log10(tabs[j]) + (random.rand() - 0.5) * 0.1);
	double eps = s_kin / 4. / E; // random background photon energy
	
	// Use approximation from A. Mastichiadis et al., Astroph. Journ. 300:178-189 (1986), eq. 30.
	// This approx is valid only for alpha >=100 where alpha = p0*eps*costheta - E0*eps
	// For our purposes, me << E0 --> p0~E0 --> alpha = E0*eps*(costheta - 1) >= 100
	double Epp = 5.7e-1 * pow(eps / mec2, -0.56) * pow(E / mec2, 0.44) * mec2;
	
	double f = Epp / E;
	
	if (haveElectrons) {
		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
		if (random.rand() < pow(1 - f, thinning)) {
			double w = 1. / pow(1 - f, thinning);
			candidate->addSecondary(11, Epp / (1 + z), pos, w, interactionTag);
		}
		if (random.rand() < pow(f, thinning)) {
			double w = 1. / pow(f, thinning);
			candidate->addSecondary(-11, Epp / (1 + z), pos, w, interactionTag);
		}
	}
	
	// Update the primary particle energy.
	// This is done after adding the secondaries to correctly set the secondaries parent
	candidate->current.setEnergy((E - 2 * Epp) / (1. + z));
	
}

void EMTripletPairProduction::process(Candidate *candidate) const {
	// check if electron / positron
	int id = candidate->current.getId();
	if (abs(id) != 11)
		return;
	
	// scale the particle energy instead of background photons
	double z = candidate->getRedshift();
	double E = (1 + z) * candidate->current.getEnergy();
	Vector3d position = candidate->current.getPosition();
	
	double scaling = pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);
	double rate = this->interactionRates->getProcessRate(E, position);
	
	if (rate < 0)
		return;
	
	rate *= scaling;
	
	// run this loop at least once to limit the step size
	double step = candidate->getCurrentStep();
	Random &random = Random::instance();
	do {
		double randDistance = -log(random.rand()) / rate;
		// check for interaction; if it doesn't occur, limit next step
		if (step < randDistance) {
			candidate->limitNextStep(limit / rate);
			return;
		}
		performInteraction(candidate);
		step -= randDistance;
	} while (step > 0.);
	
}

void EMTripletPairProduction::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string EMTripletPairProduction::getInteractionTag() const {
	return interactionTag;
}

} // namespace crpropa

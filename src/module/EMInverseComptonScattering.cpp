#include "crpropa/module/EMInverseComptonScattering.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"
#include "crpropa/Common.h"

#include "crpropa/InteractionRates.h"

#include <fstream>
#include <locale>
#include <iomanip>  // Required for std::fixed and std::setprecision
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

static const double mec2 = mass_electron * c_squared;

EMInverseComptonScattering::EMInverseComptonScattering(ref_ptr<PhotonField> photonField, bool havePhotons, double thinning, double limit, ref_ptr<Surface> surface) {

	setSurface(surface);
	setPhotonField(photonField);
	setHavePhotons(havePhotons);
	setLimit(limit);
	setThinning(thinning);
	
}

void EMInverseComptonScattering::setPhotonField(ref_ptr<PhotonField> photonField) {
	
	this->photonField = photonField;
	std::string fname = photonField->getFieldName();
	setDescription("EMInverseComptonScattering: " + fname);
	
	if (!this->photonField->hasPositionDependence()) {
		this->interactionRates = new InteractionRatesHomogeneous(
			getDataPath("EMInverseComptonScattering/rate_" + fname + ".txt"),
			getDataPath("EMInverseComptonScattering/cdf_" + fname + ".txt")
		);

	} else {
		this->interactionRates = new InteractionRatesPositionDependent(
			getDataPath("EMInverseComptonScattering/" + fname + "/Rate/"),
			getDataPath("EMInverseComptonScattering/" + fname + "/CumulativeRate/"),
			this->surface
		);
	}
}

void EMInverseComptonScattering::setHavePhotons(bool havePhotons) {
	this->havePhotons = havePhotons;
}

void EMInverseComptonScattering::setLimit(double limit) {
	this->limit = limit;
}

void EMInverseComptonScattering::setThinning(double thinning) {
	this->thinning = thinning;
}

void EMInverseComptonScattering::setSurface(ref_ptr<Surface> surface) {
		this->surface = surface;
}

ref_ptr<Surface> EMInverseComptonScattering::getSurface() const {
		return this->surface;
}

void EMInverseComptonScattering::setInteractionRates(ref_ptr<InteractionRates> intRates) {
	this->interactionRates = intRates;
}

ref_ptr<InteractionRates> EMInverseComptonScattering::getInteractionRates() const {
	return this->interactionRates;
}

void EMInverseComptonScattering::initRate(std::string path) {
	this->interactionRates->initRate(path);
}

void EMInverseComptonScattering::initCumulativeRate(std::string path) {
	this->interactionRates->initCumulativeRate(path);
}

// Class to calculate the energy distribution of the ICS photon and to sample from it
class ICSSecondariesEnergyDistribution {
	private:
		std::vector< std::vector<double> > data;
		std::vector<double> s_values;
		size_t Ns;
		size_t Nrer;
		double s_min;
		double s_max;
		double dls;

	public:
		// differential cross-section, see Lee '96 (arXiv:9604098), eq. 23 for x = Ee'/Ee
		double dSigmadE(double x, double beta) {
			double q = ((1 - beta) / beta) * (1 - 1./x);
			return ((1 + beta) / beta) * (x + 1./x + 2 * q + q * q);
		}

		// create the cumulative energy distribution of the up-scattered photon
		ICSSecondariesEnergyDistribution() {
			Ns = 1000;
			Nrer = 1000;
			s_min = mec2 * mec2;
			s_max = 2e23 * eV * eV;
			dls = (log(s_max) - log(s_min)) / Ns;
			data = std::vector< std::vector<double> >(1000, std::vector<double>(1000));
			std::vector<double> data_i(1000);

			// tabulate s bin borders
			s_values = std::vector<double>(1001);
			for (size_t i = 0; i < Ns + 1; ++i)
				s_values[i] = s_min * exp(i*dls);


			// for each s tabulate cumulative differential cross section
			for (size_t i = 0; i < Ns; i++) {
				double s = s_min * exp((i+0.5) * dls);
				double beta = (s - s_min) / (s + s_min);
				double x0 = (1 - beta) / (1 + beta);
				double dlx = -log(x0) / Nrer;

				// cumulative midpoint integration
				data_i[0] = dSigmadE(x0, beta) * expm1(dlx);
				for (size_t j = 1; j < Nrer; j++) {
					double x = x0 * exp((j+0.5) * dlx);
					double dx = exp((j+1) * dlx) - exp(j * dlx);
					data_i[j] = dSigmadE(x, beta) * dx;
					data_i[j] += data_i[j-1];
				}
				data[i] = data_i;
			}
		}

		// draw random energy for the up-scattered photon Ep(Ee, s)
		double sample(double Ee, double s) {
			size_t idx = std::lower_bound(s_values.begin(), s_values.end(), s) - s_values.begin();
			std::vector<double> s0 = data[idx];
			Random &random = Random::instance();
			size_t j = random.randBin(s0) + 1; // draw random bin (upper bin boundary returned)
			double beta = (s - s_min) / (s + s_min);
			double x0 = (1 - beta) / (1 + beta);
			double dlx = -log(x0) / Nrer;
			double binWidth = x0 * (exp(j * dlx) - exp((j-1) * dlx));
			double Ep = (x0 * exp((j-1) * dlx) + binWidth) * Ee;
			return std::min(Ee, Ep); // prevent Ep > Ee from numerical inaccuracies
		}
};

void EMInverseComptonScattering::performInteraction(Candidate *candidate) const {

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
	
	// sample the value of s
	Random &random = Random::instance();
	size_t i = closestIndex(E, tabE);
	size_t j = random.randBin(tabCDF[i]);
	double s_kin = pow(10, log10(tabs[j]) + (random.rand() - 0.5) * 0.1);
	double s = s_kin + mec2 * mec2;
	
	// sample electron energy after scattering
	static ICSSecondariesEnergyDistribution distribution;
	double Enew = distribution.sample(E, s);
	
	// add up-scattered photon
	double Esecondary = E - Enew;
	double f = Enew / E;
	if (havePhotons) {
		if (random.rand() < pow(1 - f, thinning)) {
			double w = 1. / pow(1 - f, thinning);
			Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
			candidate->addSecondary(22, Esecondary / (1 + z), pos, w, interactionTag);
		}
	}
	
	// update the primary particle energy; do this after adding the secondary to correctly set the secondary's parent
	candidate->current.setEnergy(Enew / (1 + z));
	
}

void EMInverseComptonScattering::process(Candidate *candidate) const {
	// check if electron / positron
	int id = candidate->current.getId();
	if (abs(id) != 11)
		return;
	
	// scale the particle energy instead of background photons
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);
	Vector3d position = candidate->current.getPosition();
	
	// interaction rate
	double rate = this->interactionRates->getProcessRate(E, position);
	
	if (rate < 0)
		return;
	
	rate *= pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);
	
	// run this loop at least once to limit the step size
	double step = candidate->getCurrentStep();
	Random &random = Random::instance();
	do {
		double randDistance = -log(random.rand()) / rate;
		
		// check for interaction; if it doesn't ocurr, limit next step
		if (step < randDistance) {
			candidate->limitNextStep(limit / rate);
			return;
		}
		performInteraction(candidate);
		
		// repeat with remaining step
		step -= randDistance;
	} while (step > 0);
	
}

void EMInverseComptonScattering::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string EMInverseComptonScattering::getInteractionTag() const {
	return interactionTag;
}

} // namespace crpropa

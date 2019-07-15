#include "crpropa/module/EMInverseComptonScattering.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"
#include "crpropa/Common.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

static const double mec2 = mass_electron * c_squared;

EMInverseComptonScattering::EMInverseComptonScattering(PhotonField photonField,
													   bool havePhotons,
													   std::string tag,
													   double limit) {
	setPhotonField(photonField);
	this->spaceTimeGrid = ScalarGrid4d();
	this->spaceGrid = ScalarGrid();
	this->havePhotons = havePhotons;
	this->tag = tag;
	this->limit = limit;
	setDescription("EMInverseComptonScattering_isotropicConstant");
}

EMInverseComptonScattering::EMInverseComptonScattering(PhotonField photonField,
													   ScalarGrid4d spaceTimeGrid,
													   bool havePhotons,
													   std::string tag,
													   double limit) {
	setPhotonField(photonField);
	this->spaceTimeGrid = spaceTimeGrid;
	this->spaceGrid = ScalarGrid();
	this->havePhotons = havePhotons;
	this->tag = tag;
	this->limit = limit;
	setDescription("EMInverseComptonScattering_spaceTimeDependent");
}

EMInverseComptonScattering::EMInverseComptonScattering(PhotonField photonField,
													   ScalarGrid spaceGrid,
													   bool havePhotons,
													   std::string tag,
													   double limit) {
	setPhotonField(photonField);
	this->spaceTimeGrid = ScalarGrid4d();
	this->spaceGrid = spaceGrid;
	this->havePhotons = havePhotons;
	this->tag = tag;
	this->limit = limit;
	setDescription("EMInverseComptonScattering_spaceDependentConstant");
}

void EMInverseComptonScattering::setPhotonField(PhotonField photonField) {
	this->photonField = photonField;
	std::string fname = photonFieldName(photonField);
	initRate(getDataPath("EMInverseComptonScattering/rate_" + fname + ".txt"));
	initCumulativeRate(getDataPath("EMInverseComptonScattering/cdf_" + fname + ".txt"));
}

void EMInverseComptonScattering::setHavePhotons(bool havePhotons) {
	this->havePhotons = havePhotons;
}

void EMInverseComptonScattering::setLimit(double limit) {
	this->limit = limit;
}

void EMInverseComptonScattering::initRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error("EMInverseComptonScattering: could not open file " + filename);

	// clear previously loaded tables
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

void EMInverseComptonScattering::initCumulativeRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error("EMInverseComptonScattering: could not open file " + filename);

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
			s_max = 1e23 * eV * eV;
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
				data_i[0] = dSigmadE(x0, beta) * (exp(dlx) - 1);
				for (size_t j = 1; j < Nrer; j++) {
					double x = x0 * exp((j+0.5)*dlx);
					double dx = exp((j+1)*dlx) - exp(j*dlx);
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
			double binWidth = x0 * (exp(j*dlx) - exp((j-1)*dlx));
			double Ep = (x0 * exp((j-1)*dlx) + binWidth) * Ee;
			return std::min(Ee, Ep); // prevent Ep > Ee from numerical inaccuracies
		}
};

void EMInverseComptonScattering::performInteraction(Candidate *candidate) const {
	// scale the particle energy instead of background photons
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	if (E < tabE.front() or E > tabE.back())
		return;

	// geometric scaling
	Vector3d pos = candidate->current.getPosition();
    const double time = candidate->getTrajectoryLength() / c_light;

	double geometricScaling = 1.;
	const std::string description = getDescription();
	if (description == "EMInverseComptonScattering_isotropicConstant") {
		// do nothing, just check for correct initialization
	} else if (description == "EMInverseComptonScattering_spaceDependentConstant") {
		geometricScaling *= spaceGrid.interpolate(pos);
	} else if (description == "EMInverseComptonScattering_spaceTimeDependent") {
		geometricScaling *= spaceTimeGrid.interpolate(pos, time);
	} else {
		throw std::runtime_error("EMInverseComptonScattering: invalid description string");
	}
	if (geometricScaling == 0.)
		return;

	// sample the value of s
	Random &random = Random::instance();
	size_t i = closestIndex(E, tabE);
	size_t j = random.randBin(tabCDF[i]);
	double s_kin = pow(10, log10(tabs[j] * geometricScaling) + (random.rand() - 0.5) * 0.1);
	double s = s_kin + mec2 * mec2;

	// sample electron energy after scattering
	static ICSSecondariesEnergyDistribution distribution;
	double Enew = distribution.sample(E, s);

	// add up-scattered photon
	double Esecondary = E - Enew;
	if (havePhotons) {
		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
		candidate->addSecondary(22, Esecondary / (1 + z), pos, tag);
	}

	// update the primary particle energy; do this after adding the secondary to correctly set the secondary's parent
	candidate->current.setEnergy(Enew / (1 + z));
}

void EMInverseComptonScattering::process(Candidate *candidate) const {
	// check if electron / positron
	int id = candidate->current.getId();
	if (id != 11 && id != -11)
		return;

	// scale the particle energy instead of background photons
	double z = candidate->getRedshift();
	double E = (1 + z) * candidate->current.getEnergy();

	if (E < tabEnergy.front() or (E > tabEnergy.back()))
		return;

	// interaction rate
	// geometric scaling
	Vector3d pos = candidate->current.getPosition();
    const double time = candidate->getTrajectoryLength() / c_light;

    double rate = 1.;
	const std::string description = getDescription();
	if (description == "EMInverseComptonScattering_isotropicConstant") {
		// do nothing, just check for correct initialization
	} else if (description == "EMInverseComptonScattering_spaceDependentConstant") {
		rate *= spaceGrid.interpolate(pos);
	} else if (description == "EMInverseComptonScattering_spaceTimeDependent") {
		rate *= spaceTimeGrid.interpolate(pos, time);
	} else {
		throw std::runtime_error("EMInverseComptonScattering: invalid description string");
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

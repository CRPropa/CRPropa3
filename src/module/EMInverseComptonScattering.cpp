#include "crpropa/module/EMInverseComptonScattering.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

EMInverseComptonScattering::EMInverseComptonScattering(PhotonField photonField,
		bool havePhotons, double limit) {
	setPhotonField(photonField);
	this->havePhotons = havePhotons;
	this->limit = limit;
}

void EMInverseComptonScattering::setPhotonField(PhotonField photonField) {
	this->photonField = photonField;
	switch (photonField) {
	case CMB:
		setDescription("EMInverseComptonScattering: CMB");
		initRate(getDataPath("EMInverseComptonScattering_CMB.txt"));
		initCumulativeRate(getDataPath("EMInverseComptonScattering_CDF_CMB.txt"));
		break;
	case IRB:  // default: Kneiske '04 IRB model
	case IRB_Kneiske04:
		setDescription("EMInverseComptonScattering: IRB (Kneiske 2004)");
		initRate(getDataPath("EMInverseComptonScattering_IRB_Kneiske04.txt"));
		initCumulativeRate(getDataPath("EMInverseComptonScattering_CDF_IRB_Kneiske04.txt"));
		break;
	case IRB_Stecker05:
		setDescription("EMInverseComptonScattering: IRB (Stecker 2005)");
		initRate(getDataPath("EMInverseComptonScattering_IRB_Stecker05.txt"));
		initCumulativeRate(getDataPath("EMInverseComptonScattering_CDF_IRB_Stecker05.txt"));
		break;
	case IRB_Franceschini08:
		setDescription("EMInverseComptonScattering: IRB (Franceschini 2008)");
		initRate(getDataPath("EMInverseComptonScattering_IRB_Franceschini08.txt"));
		initCumulativeRate(getDataPath("EMInverseComptonScattering_CDF_IRB_Franceschini08.txt"));
		break;
	case IRB_Finke10:
		setDescription("EMInverseComptonScattering: IRB (Finke 2010)");
		initRate(getDataPath("EMInverseComptonScattering_IRB_Finke10.txt"));
		initCumulativeRate(getDataPath("EMInverseComptonScattering_CDF_IRB_Finke10.txt"));
		break;
	case IRB_Dominguez11:
		setDescription("EMInverseComptonScattering: IRB (Dominguez 2011)");
		initRate(getDataPath("EMInverseComptonScattering_IRB_Dominguez11.txt"));
		initCumulativeRate(getDataPath("EMInverseComptonScattering_CDF_IRB_Dominguez11.txt"));
		break;
	case IRB_Gilmore12:
		setDescription("EMInverseComptonScattering: IRB (Gilmore 2012)");
		initRate(getDataPath("EMInverseComptonScattering_IRB_Gilmore12.txt"));
		initCumulativeRate(getDataPath("EMInverseComptonScattering_CDF_IRB_Gilmore12.txt"));
		break;
	case URB_Protheroe96:
		setDescription("EMInverseComptonScattering: URB (Protheroe 1996)");
		initRate(getDataPath("EMInverseComptonScattering_URB_Protheroe96.txt"));
		initCumulativeRate(getDataPath("EMInverseComptonScattering_CDF_URB_Protheroe96.txt"));
		break;
	default:
		throw std::runtime_error(
				"EMInverseComptonScattering: unknown photon background");
	}
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
		throw std::runtime_error(
				"EMInverseComptonScattering: could not open file " + filename);

	// clear previously loaded interaction rates
	tabElectronEnergy.clear();
	tabInteractionRate.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabElectronEnergy.push_back(pow(10, a) * eV);
				tabInteractionRate.push_back(b / Mpc);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
}

void EMInverseComptonScattering::initCumulativeRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error(
				"EMInverseComptonScattering: could not open file " + filename);

	// clear previously loaded interaction rates
	tabE.clear();
	tabs.clear();
	tabCumulativeRate.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b, c;
			infile >> a >> b >> c;
			if (infile) {
				tabE.push_back(pow(10, a) * eV);
				tabs.push_back(pow(10,b) * eV * eV);
				tabCumulativeRate.push_back(c / Mpc);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
}

// Differential cross-section for inverse Compton scattering. from Lee, eq. 23 for x = Ee'/Ee
// compare Lee 96 arXiv:9604098
double dSigmadE_ICSx(double x, double beta) {
	double q = ((1 - beta) / beta) * (1 - 1./x);
	double A = x + 1./x;
	return ((1 + beta) / beta) * (A + 2 * q + q * q);
}

// Hold an data array to interpolate the energy distribution on 
class ICSSecondariesEnergyDistribution {
	private:
		double ElectronMass;
		std::vector< std::vector<double> > data;
		std::vector<double> s_values;
		size_t Ns;
		size_t Nrer;
		double s_min;
		double s_max;
		double dls;

	public:
		ICSSecondariesEnergyDistribution() {
			ElectronMass = mass_electron * c_squared;
			Ns = 1000;
			Nrer = 1000;
			s_min = ElectronMass * ElectronMass;
			s_max = 1e23 * eV * eV;
			dls = (log(s_max) - log(s_min)) / Ns;
			data = std::vector< std::vector<double> >(1000,std::vector<double>(1000));
			s_values = std::vector<double>(1001);

			for (size_t i = 0; i < Ns + 1; ++i)
				s_values[i] = s_min * exp(i*dls); // tabulate s bin borders

			for (size_t i = 0; i < Ns; i++) {
				const double s = s_min * exp(i*dls + 0.5*dls); // choose bin centers for evaluation to stay away from critical lower boundary s = m_e**2
				double beta = (s - ElectronMass * ElectronMass) / (s + ElectronMass * ElectronMass);
				double x0 = (1.-beta) / (1.+beta);
				double dx = -log((1. - beta)/(1.+beta)) / Nrer;
				std::vector<double> data_i(1000);
				data_i[0] = dSigmadE_ICSx(x0,beta)*(exp(dx)-1.);
				for (size_t j = 1; j < Nrer; j++){
					double x = x0 * exp(j*dx + 0.5*dx); 
					data_i[j] = dSigmadE_ICSx(x, beta)*(exp((j+1)*dx)-exp(j*dx)) + data_i[j-1]; // cumulative midpoint integration 
				}
				data[i] = data_i;
			}
		}

		// returns pointer to the the integrated distribution for a given s
		std::vector<double> getDistribution(double s) {
			size_t idx = std::lower_bound(s_values.begin(), s_values.end(), s) - s_values.begin();
			std::vector<double> s0 = data[idx];
			return s0;
		}

		// samples the integrated distribution and returns Eer(Ee, s)
		double sample(double Ee, double s) {
			std::vector<double> s0 = getDistribution(s); 
			Random &random = Random::instance();
			size_t j = random.randBin(s0) + 1; // draw random bin (upper bin boundary returned), cause 0 not in CDF +1 to get right x value
			double beta = (s - ElectronMass * ElectronMass) / (s + ElectronMass * ElectronMass);
			double x0 = (1-beta) / (1+beta);
			double dx =  - log((1-beta) / (1+beta)) / Nrer;
			double binWidth = x0*(exp(j*dx)-exp((j-1)*dx));
			return (x0*exp((j-1)*dx) + binWidth) * Ee; 
		}
};

void EMInverseComptonScattering::performInteraction(Candidate *candidate) const {

	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy();
	double mec2 = mass_electron * c_squared;

	// interpolate between tabulated electron energies to get corresponding cdf
	if (E > tabE.back() || E < tabE.front())
		return;
	size_t i = std::upper_bound(tabE.begin(), tabE.end(), E) - tabE.begin() - 500;
	double a = (E - tabE[i]) / (tabE[i + 500] - tabE[i]);

	std::vector<double> cdf(500);
	for (size_t j = 0; j < 500; j++)
		cdf[j] = tabCumulativeRate[i+j] + a * (tabCumulativeRate[i+500+j] - tabCumulativeRate[i+j]);

	// draw random value between 0. and maximum of corresponding cdf
	// choose bin of s where cdf(s) = cdf_rand -> s_rand
	Random &random = Random::instance();
	size_t j = random.randBin(cdf); // draw random bin (upper bin boundary returned)
	double binWidth = (tabs[i+j] - tabs[i+j-1]);
	double s_kin = tabs[i+j-1] + random.rand() * binWidth; // draw random s uniformly distributed in bin
	s_kin *= (1 + z);
	double s = s_kin + mec2*mec2;

	static ICSSecondariesEnergyDistribution interpolation;
	double Epost = interpolation.sample(E,s);
	if (havePhotons) {
		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(),candidate->current.getPosition());
		candidate->addSecondary(22, (E-Epost), pos);
	}
	candidate->current.setEnergy(Epost);
}

void EMInverseComptonScattering::process(Candidate *candidate) const {
	double step = candidate->getCurrentStep();
	double z = candidate->getRedshift();

	// check if electron / positron
	int id = candidate->current.getId();
	if (id != 11 && id != -11)
		return;

	// instead of scaling the background photon energies, scale the electron energy
	double E = (1 + z) * candidate->current.getEnergy();

	// check if in tabulated energy range
	if (E < tabElectronEnergy.front() or (E > tabElectronEnergy.back()))
		return;

	// find interaction with minimum random distance
	Random &random = Random::instance();
	double randDistance = std::numeric_limits<double>::max();

	// cosmological scaling of interaction distance (comoving)
	double scaling = pow(1 + z, 2) * photonFieldScaling(photonField, z);
	double rate = scaling * interpolate(E, tabElectronEnergy, tabInteractionRate);
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

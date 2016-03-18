#include "crpropa/module/EMPairProduction.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <stdexcept>


namespace crpropa {

EMPairProduction::EMPairProduction(PhotonField photonField,
		bool haveElectrons, double limit) {
	setPhotonField(photonField);
	this->haveElectrons = haveElectrons;
	this->limit = limit;
}

void EMPairProduction::setPhotonField(PhotonField photonField) {
	this->photonField = photonField;
	switch (photonField) {
	case CMB:
		setDescription("EMPairProduction: CMB");
		initRate(getDataPath("EMPairProduction_CMB.txt"));
		initCumulativeRate(getDataPath("EMPairProduction_CDF_CMB.txt"));
		break;
	case IRB:  // default: Kneiske '04 IRB model
	case IRB_Kneiske04:
		setDescription("EMPairProduction: IRB (Kneiske 2004)");
		initRate(getDataPath("EMPairProduction_IRB_Kneiske04.txt"));
		initCumulativeRate(getDataPath("EMPairProduction_CDF_IRB_Kneiske04.txt"));
		break;
	case IRB_Stecker05:
		setDescription("EMPairProduction: IRB (Stecker 2005)");
		initRate(getDataPath("EMPairProduction_IRB_Stecker05.txt"));
		initCumulativeRate(getDataPath("EMPairProduction_CDF_IRB_Stecker05.txt"));
		break;
	case IRB_Franceschini08:
		setDescription("EMPairProduction: IRB (Franceschini 2008)");
		initRate(getDataPath("EMPairProduction_IRB_Franceschini08.txt"));
		initCumulativeRate(getDataPath("EMPairProduction_CDF_IRB_Franceschini08.txt"));
		break;
	case IRB_Finke10:
		setDescription("EMPairProduction: IRB (Finke 2010)");
		initRate(getDataPath("EMPairProduction_IRB_Finke10.txt"));
		initCumulativeRate(getDataPath("EMPairProduction_CDF_IRB_Finke10.txt"));
		break;
	case IRB_Dominguez11:
		setDescription("EMPairProduction: IRB (Dominguez 2011)");
		initRate(getDataPath("EMPairProduction_IRB_Dominguez11.txt"));
		initCumulativeRate(getDataPath("EMPairProduction_CDF_IRB_Dominguez11.txt"));
		break;
	case IRB_Gilmore12:
		setDescription("EMPairProduction: IRB (Gilmore 2012)");
		initRate(getDataPath("EMPairProduction_IRB_Gilmore12.txt"));
		initCumulativeRate(getDataPath("EMPairProduction_CDF_IRB_Gilmore12.txt"));
		break;
	case URB_Protheroe96:
		setDescription("EMPairProduction: URB (Protheroe 1996)");
		initRate(getDataPath("EMPairProduction_URB_Protheroe96.txt"));
		initCumulativeRate(getDataPath("EMPairProduction_CDF_URB_Protheroe96.txt"));
		break;
	default:
		throw std::runtime_error(
				"EMPairProduction: unknown photon background");
	}
}

void EMPairProduction::setHaveElectrons(bool haveElectrons) {
	this->haveElectrons = haveElectrons;
}

void EMPairProduction::setLimit(double limit) {
	this->limit = limit;
}

void EMPairProduction::initRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error(
				"EMPairProduction: could not open file " + filename);

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

void EMPairProduction::initCumulativeRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error(
				"EMPairProduction: could not open file " + filename);

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
				tabs.push_back(pow(10,b) * eV*eV);
				tabCumulativeRate.push_back(c / Mpc);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
}

///	 Differential cross-section for pair production for x = Epositron/Egamma
double dSigmadE_PPx(double x, double beta) {
	const double A = (x / (1. - x) + (1. - x) / x );
	const double B =  (1. / x + 1. / (1. - x) );
	return A + (1. - beta*beta) * B - (1. - beta*beta) * (1. - beta*beta) / 4. * B*B;
}

/// Hold an data array to interpolate the energy distribution on
class PPSecondariesEnergyDistribution
{
	private:
		double *_data;
		size_t _Ns;
		size_t _Nrer;
		double _s_min;
		double _s_max;
		double _dls;

	public:
		PPSecondariesEnergyDistribution(double s_min = 4. * mass_electron*c_squared * mass_electron*c_squared, double s_max =1e23,
				size_t Ns = 1000, size_t Nrer = 1000 )
		{
			if (s_min < 4.*mass_electron*c_squared * mass_electron*c_squared)
			{
				std::cerr << "Warning: Minimum COM Energy in PP Interpolation s = " << s_min << " <  (2*m_e)**2 selected. Setting to s_min = (2*m_e)**2.\n";
				s_min = 4.*mass_electron*c_squared * mass_electron*c_squared;
			}
			_Ns = Ns;
			_Nrer = Nrer;
			_s_min =s_min;
			_s_max = s_max;
			_data = new double[Ns*Nrer];
			_dls = (log(s_max) - log(s_min)) / (Ns);
			double ElectronMass = mass_electron*c_squared;

			for (size_t i = 0; i < Ns; i++)
			{
				const double s = s_min * exp(i*_dls);
				double beta = sqrt(1. - 4. * ElectronMass*ElectronMass /s);
				
				double x0 = log((1.-beta) / 2.);
				double dx = ( log((1. + beta)/2) -  log((1.-beta) / 2.)) / (Nrer); 
				_data[i *Nrer] = 0;
				for (size_t j = 1; j < Nrer; j++)
				{
					double x = exp(x0 + j*dx); 
					_data[i * Nrer + j] =	dSigmadE_PPx(x, beta) + _data[i*Nrer +j-1]; //TODO: tables contains some nans
				}
			}
		}

		// returns pointer to the the integrated distribution for a given s
		double* getDistribution(double s)
		{
			size_t idx = (log(s / _s_min)) / _dls;
			double *s0 = &_data[idx * _Nrer];
			return s0;
		}

		//samples the integrated distribution and returns Eer(Ee, s)
		double sample(double E0,double s)
		{
			double ElectronMass = mass_electron*c_squared;
			double *s0 = getDistribution(s); 
			Random &random = Random::instance();
			double rnd = random.rand() *s0[_Nrer-1];

			for (size_t i=0; i < _Nrer; i++)
			{
				if (rnd < s0[i])
				{
					double beta = sqrt(1. - 4.* ElectronMass * ElectronMass / s);

					double x0 = log((1.-beta) / 2.);
					double dx = ( log((1. + beta)/2) -  log((1.-beta) / 2.)) / (_Nrer);
					if (random.rand() < 0.5)
						return exp(x0 + (i)*dx) * E0;
					else
						return E0 * (1-exp(x0 + (i)*dx));
				}
			}
			std::cerr << "PPSecondariesEnergyDistribution out of bounds!" << std::endl;
			std::cerr << "  s0[0] = " << s0[0] << "  s0[_Nrer-1] = " << s0[_Nrer-1] << "  rnd = " << rnd << std::endl;
			throw std::runtime_error("Grave logic error in PPSecondariesEnergyDistribution!");
		}
};

// Helper function for actual Monte Carlo sampling to avoid code-duplication
double __extractPPSecondariesEnergy(double E0, double s)
{
	static PPSecondariesEnergyDistribution interpolation;
	return interpolation.sample(E0, s);
}

void EMPairProduction::performInteraction(Candidate *candidate) const {

	if (haveElectrons){
		int id = candidate->current.getId();
		double z = candidate->getRedshift();
		double E = candidate->current.getEnergy();
		double Epos = 0.;
		double mec2 = mass_electron * c_squared;

		// interpolate between tabulated electron energies to get corresponding cdf
		size_t i = std::upper_bound(tabE.begin(), tabE.end(), E) - tabE.begin() - 500;
		double a = (E - tabE[i]) / (tabE[i + 500] - tabE[i]);
		if (E < tabE.front() || E > tabE.back())
			return;

		std::vector<double> cdf(500);
		for (size_t j = 0; j < 500; j++)
			cdf[j] = tabCumulativeRate[i+j] + a * (tabCumulativeRate[i+500+j] - tabCumulativeRate[i+j]);

		// draw random value between 0. and maximum of corresponding cdf
		// choose bin of s where cdf(s) = cdf_rand -> s_rand
		Random &random = Random::instance();
		size_t j = random.randBin(cdf); // draw random bin
		double binWidth = (tabs[i+j+1] - tabs[i+j]);
		double s_kin = tabs[i+j] + random.rand() * binWidth; // draw random s uniformly distributed in bin
		s_kin *(1 + z);
		if (s_kin < 4*mec2*mec2)
			std::cout << "ERROR" << std::endl;
		Epos = __extractPPSecondariesEnergy(E,s_kin);

		Vector3d pos = randomPositionInPropagationStep(candidate);
		candidate->addSecondary(-11, (E-Epos), pos);
		candidate->addSecondary(11, Epos, pos);
	}
	candidate->setActive(false);
}

void EMPairProduction::process(Candidate *candidate) const {
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

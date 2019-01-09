#include "crpropa/module/EMPairProduction.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <stdexcept>


namespace crpropa {

static const double mec2 = mass_electron * c_squared;

EMPairProduction::EMPairProduction(PhotonField photonField,
								   ScalarGrid4d geometryGrid,
								   bool haveElectrons,
								   double limit) : haveElectrons(haveElectrons), 
												   limit(limit) {
	setPhotonField(photonField);
	this->geometryGrid = geometryGrid;
}

void EMPairProduction::setPhotonField(PhotonField photonField) {
	this->photonField = photonField;
	std::string fname = photonFieldName(photonField);
	setDescription("EMPairProduction: " + fname);
	initRate(getDataPath("EMPairProduction/rate_" + fname + ".txt"));
	initCumulativeRate(getDataPath("EMPairProduction/cdf_" + fname + ".txt"));
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
		throw std::runtime_error("EMPairProduction: could not open file " + filename);

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

void EMPairProduction::initCumulativeRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error("EMPairProduction: could not open file " + filename);

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

// Hold an data array to interpolate the energy distribution on
class PPSecondariesEnergyDistribution {
	private:
		std::vector<double> tab_s;
		std::vector< std::vector<double> > data;
		size_t N;

	public:
		// differential cross section for pair production for x = Epositron/Egamma, compare Lee 96 arXiv:9604098
		double dSigmadE_PPx(double x, double beta) {
			double A = (x / (1. - x) + (1. - x) / x );
			double B =  (1. / x + 1. / (1. - x) );
			double y = (1 - beta * beta);
			return A + y * B - y * y / 4 * B * B;
		}

		PPSecondariesEnergyDistribution() {
			N = 1000;
			size_t Ns = 1000;
			double s_min = 4 * mec2 * mec2;
			double s_max = 1e23 * eV * eV;
			double dls = log(s_max / s_min) / Ns;
			data = std::vector< std::vector<double> >(Ns, std::vector<double>(N));
			tab_s = std::vector<double>(Ns + 1);

			for (size_t i = 0; i < Ns + 1; ++i)
				tab_s[i] = s_min * exp(i*dls); // tabulate s bin borders

			for (size_t i = 0; i < Ns; i++) {
				double s = s_min * exp(i*dls + 0.5*dls);
				double beta = sqrt(1 - s_min/s);
				double x0 = (1 - beta) / 2;
				double dx = log((1 + beta) / (1 - beta)) / N;

				// cumulative midpoint integration
				std::vector<double> data_i(1000);
				data_i[0] = dSigmadE_PPx(x0, beta) * (exp(dx) - 1);
				for (size_t j = 1; j < N; j++) {
					double x = x0 * exp(j*dx + 0.5*dx);
					double binWidth = exp((j+1)*dx)-exp(j*dx);
					data_i[j] = dSigmadE_PPx(x, beta) * binWidth + data_i[j-1];
				}
				data[i] = data_i;
			}
		}

		// sample positron energy from cdf(E, s_kin)
		double sample(double E0, double s) {
			// get distribution for given s
			size_t idx = std::lower_bound(tab_s.begin(), tab_s.end(), s) - tab_s.begin();
			std::vector<double> s0 = data[idx];

			// draw random bin
			Random &random = Random::instance();
			size_t j = random.randBin(s0) + 1;

			double s_min = 4 * mec2 * mec2;
			double beta = sqrt(1 - s_min / s);
			double x0 = (1 - beta) / 2.;
			double dx = log((1 + beta) / (1 - beta)) / N;
			double binWidth = x0 * (exp(j*dx) - exp((j-1)*dx));
			if (random.rand() < 0.5)
				return E0 * (x0 * exp((j-1) * dx) + binWidth);
			else
				return E0 * (1 - (x0 * exp((j-1) * dx) + binWidth));
		}
};

void EMPairProduction::performInteraction(Candidate *candidate) const {
	// scale particle energy instead of background photon energy
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	// cosmic ray photon is lost after interacting
	candidate->setActive(false);

	// check if secondary electron pair needs to be produced
	if (not haveElectrons)
		return;

	// check if in tabulated energy range
	if (E < tabE.front() or (E > tabE.back()))
		return;

	// sample the value of s
	Random &random = Random::instance();
	size_t i = closestIndex(E, tabE);  // find closest tabulation point

	// geometric scaling
	Vector3d pos = candidate->current.getPosition();
    double time = candidate->getTrajectoryLength()/c_light;
    double geometricScaling = geometryGrid.interpolate(pos, time);
    if (geometricScaling == 0.)
    	return;

    std::vector<double> tabCDF_geoScaled;
    for (int j = 0; j < tabCDF[i].size(); ++j) {
    	tabCDF_geoScaled.push_back(tabCDF[i][j]*geometricScaling);
    }
	size_t j = random.randBin(tabCDF_geoScaled);
	double lo = std::max(4 * mec2 * mec2, tabs[j-1]);  // first s-tabulation point below min(s_kin) = (2 me c^2)^2; ensure physical value
	double hi = tabs[j];
	double s = lo + random.rand() * (hi - lo);

	// sample electron / positron energy
	static PPSecondariesEnergyDistribution interpolation;
	double Ee = interpolation.sample(E, s);
	double Ep = E - Ee;

	// sample random position along current step
	pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
	candidate->addSecondary(-11, Ee / (1 + z), pos);
	candidate->addSecondary(11, Ep / (1 + z), pos);
}

void EMPairProduction::process(Candidate *candidate) const {
	// check if photon
	if (candidate->current.getId() != 22)
		return;

	// scale particle energy instead of background photon energy
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	// check if in tabulated energy range
	if (E < tabEnergy.front() or (E > tabEnergy.back()))
		return;

	// interaction rate
	// geometric scaling
	Vector3d pos = candidate->current.getPosition();
    double time = candidate->getTrajectoryLength()/c_light;
	double rate = geometryGrid.interpolate(pos, time);
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

#include "crpropa/PhotonBackground.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <fstream>
#include <stdexcept>
#include <limits>
#include <cmath>

namespace crpropa {

TabularPhotonField::TabularPhotonField(std::string fieldName, bool isRedshiftDependent) {
	this->fieldName = fieldName;
	this->isRedshiftDependent = isRedshiftDependent;

	readPhotonEnergy(getDataPath("") + "Scaling/" + this->fieldName + "_photonEnergy.txt");
	readPhotonDensity(getDataPath("") + "Scaling/" + this->fieldName + "_photonDensity.txt");
	if (this->isRedshiftDependent)
		readRedshift(getDataPath("") + "Scaling/" + this->fieldName + "_redshift.txt");

	checkInputData();

	if (this->isRedshiftDependent)
		initRedshiftScaling();
}

double TabularPhotonField::getPhotonDensity(double Ephoton, double z) const {
	if (this->isRedshiftDependent) {
		return interpolate2d(Ephoton, z, this->photonEnergies, this->redshifts, this->photonDensity);
	} else {
		return interpolate(Ephoton, this->photonEnergies, this->photonDensity);
	}
}


double TabularPhotonField::getRedshiftScaling(double z) const {
	if (this->isRedshiftDependent) {
		if (z > this->redshifts.back()) {
			return 0.;
		} else if (z < this->redshifts.front()) {
			return 1.;
		} else {
			return interpolate(z, this->redshifts, this->redshiftScalings);
		}
	} else {
		return 1.;
	}
}

double TabularPhotonField::getMinimumPhotonEnergy(double z) const{
	return photonEnergies[0];
}

double TabularPhotonField::getMaximumPhotonEnergy(double z) const{
	return photonEnergies[photonEnergies.size() -1];
}

void TabularPhotonField::readPhotonEnergy(std::string filePath) {
	std::ifstream infile(filePath.c_str());
	if (!infile.good())
		throw std::runtime_error("TabularPhotonField::readPhotonEnergy: could not open " + filePath);

	std::string line;
	while (std::getline(infile, line)) {
		if (line.size() > 0)
			this->photonEnergies.push_back(std::stod(line)/eV); //conversion from [J] to [eV]
	}
	infile.close();
}

void TabularPhotonField::readPhotonDensity(std::string filePath) {
	std::ifstream infile(filePath.c_str());
	if (!infile.good())
		throw std::runtime_error("TabularPhotonField::readPhotonDensity: could not open " + filePath);

	std::string line;
	while (std::getline(infile, line)) {
		if (line.size() > 0)
			this->photonDensity.push_back(std::stod(line));
	}
	infile.close();
}

void TabularPhotonField::readRedshift(std::string filePath) {
	std::ifstream infile(filePath.c_str());
	if (!infile.good())
		throw std::runtime_error("TabularPhotonField::initRedshift: could not open " + filePath);

	std::string line;
	while (std::getline(infile, line)) {
		if (line.size() > 0)
			this->redshifts.push_back(std::stod(line));
	}
	infile.close();
}

void TabularPhotonField::initRedshiftScaling() {
	double n0 = 0.;
	for (int i = 0; i < this->redshifts.size(); ++i) {
		double z = this->redshifts[i];
		double n = 0.;
		for (int j = 0; j < this->photonEnergies.size(); ++j) {
			double e = this->photonEnergies[j];
			if (z == 0.)
				n0 += getPhotonDensity(e, z);
			n += getPhotonDensity(e, z);
		}
		this->redshiftScalings.push_back(n / n0);
	}
}

void TabularPhotonField::checkInputData() const {
	if (this->isRedshiftDependent) {
		if (this->photonDensity.size() != this->photonEnergies.size() * this-> redshifts.size())
			throw std::runtime_error("TabularPhotonField::checkInputData: length of photon density input is unequal to length of photon energy input times length of redshift input");
	} else {
		if (this->photonEnergies.size() != this->photonDensity.size())
			throw std::runtime_error("TabularPhotonField::checkInputData: length of photon energy input is unequal to length of photon density input");
	}

	for (int i = 0; i < this->photonEnergies.size(); ++i) {
		double ePrevious = 0.;
		double e = this->photonEnergies[i];
		if (e <= 0.)
			throw std::runtime_error("TabularPhotonField::checkInputData: a value in the photon energy input is not positive");
		if (e <= ePrevious)
			throw std::runtime_error("TabularPhotonField::checkInputData: photon energy values are not strictly increasing");
		ePrevious = e;
	}

	for (int i = 0; i < this->photonDensity.size(); ++i) {
		if (this->photonDensity[i] < 0.)
			throw std::runtime_error("TabularPhotonField::checkInputData: a value in the photon density input is negative");
	}

	if (this->isRedshiftDependent) {
		if (this->redshifts[0] != 0.)
			throw std::runtime_error("TabularPhotonField::checkInputData: redshift input must start with zero");

		for (int i = 0; i < this->redshifts.size(); ++i) {
			double zPrevious = -1.;
			double z = this->redshifts[i];
			if (z < 0.)
				throw std::runtime_error("TabularPhotonField::checkInputData: a value in the redshift input is negative");
			if (z <= zPrevious)
				throw std::runtime_error("TabularPhotonField::checkInputData: redshift values are not strictly increasing");
			zPrevious = z;
		}

		for (int i = 0; i < this->redshiftScalings.size(); ++i) {
			double scalingFactor = this->redshiftScalings[i];
			if (scalingFactor <= 0.)
				throw std::runtime_error("TabularPhotonField::checkInputData: initRedshiftScaling has created a non-positive scaling factor");
		}
	}
}

BlackbodyPhotonField::BlackbodyPhotonField(std::string fieldName, double blackbodyTemperature) {
	this->fieldName = fieldName;
	this->blackbodyTemperature = blackbodyTemperature;
	this->quantile = 0.0001;
}

double BlackbodyPhotonField::getPhotonDensity(double Ephoton, double z) const {
	return 8 * M_PI * pow_integer<3>(Ephoton / (h_planck * c_light)) / std::expm1(Ephoton / (k_boltzmann * this->blackbodyTemperature));
}

double BlackbodyPhotonField::getMinimumPhotonEnergy(double z) const {
	double A;
	int quantile_int = 10000*quantile;
	switch (quantile_int)
	{
	case 1:	// 0.01 % percentil
		A = 1.093586e-5 * eV / kelvin;
		break;
	case 10:		// 0.1 % percentil
		A = 2.402189e-5 * eV / kelvin; 
		break;
	case 100:		// 1 % percentil
		A = 5.417942e-5 * eV / kelvin;
		break;
	default:
		throw std::runtime_error("Quantile not understood. Please use 0.01 (1%), 0.001 (0.1%) or 0.0001 (0.01%) \n");
		break;
	}
	return A * this->blackbodyTemperature;
}

double BlackbodyPhotonField::getMaximumPhotonEnergy(double z) const{
	double A;
	int quantile_int = 10000*quantile;
	switch (quantile_int)
	{
	case 1:	// 99.99 % percentil
		A = 1.36246e-3 * eV/kelvin ;
		break;
	case 10:		// 99.9 % percentil
		A = 1.116975e-3 * eV/kelvin; 
		break;
	case 100:		// 99 % percentil
		A = 8.5562445e-4 * eV/kelvin; 
		break;
	default:
		throw std::runtime_error("Quantile not understood. Please use 0.01 (1%), 0.001 (0.1%) or 0.0001 (0.01%) \n");
		break;
	}
	return A*this->blackbodyTemperature;
}

void BlackbodyPhotonField::setQuantile(double q){
	if(not ((q==0.0001) or (q==0.001) or (q== 0.01)))
		throw std::runtime_error("Quantile not understood. Please use 0.01 (1%), 0.001 (0.1%) or 0.0001 (0.01%) \n");
	this->quantile = q;
}

PhotonFieldSampling::PhotonFieldSampling() {
	const std::string photonFieldName = "CMB";
	photonField = new CMB();
}

PhotonFieldSampling::PhotonFieldSampling(ref_ptr<PhotonField> field) {
	const std::string photonFieldName = field->getFieldName();
	photonField = field;
}

double PhotonFieldSampling::sample_eps(bool onProton, double Ein, double z) const {
	double eps = 0.;
	double epsMin = photonField->getMinimumPhotonEnergy(0)/eV;
	double epsMax = photonField->getMaximumPhotonEnergy(0)/eV;

	double pEpsMax = prob_eps_max(onProton, Ein, z, epsMin, epsMax);
	
	std::cout << "pEpsMax is done and result to: " << pEpsMax << "\n";
	// sample eps between epsMin ... epsMax
	double pEps = 0.;
	Random &random = Random::instance();
	int iRep = 0;
	const int iMax = 100000;
	do {
		iRep++;
		eps = epsMin + random.rand() * (epsMax - epsMin);
		pEps = prob_eps(eps, onProton, Ein, z);
		if (iRep > iMax)
			throw std::runtime_error("error: no photon found in sample_eps, please make sure that photon field provides photons for the interaction by adapting the energy range of the tabulated photon field.");
			break;
	} while (random.rand() * pEpsMax > pEps);
	
	return eps * eV;
}

double PhotonFieldSampling::prob_eps_max(bool onProton, double Ein, double z, double epsMin, double epsMax) const {
	// find pMax for current interaction to have a reference for the MC rejection
	double pEpsMaxTested = 0.;
	double epsDummy = 0.;
	// resolution of sampling the range between epsMin ... epsMax for finding the phtoton energy with the maximal interaction prop. 
	const int nrIteration = 100;
	for (int i = 0; i < nrIteration; ++i) {
		epsDummy = epsMin + (epsMax - epsMin) / nrIteration * i;
		const double pEpsDummy = this->prob_eps(epsDummy, onProton, Ein, z);
		if (pEpsDummy > pEpsMaxTested)
			pEpsMaxTested = pEpsDummy;
	}
	// the following factor corrects for only trying to find the maximum on nrIteration phtoton energies
	// the factor should be determined in convergence tests
	const double maxCorrectionFactor = 1.6;
	double pEpsMax = pEpsMaxTested * maxCorrectionFactor;
	return pEpsMax;
}

double PhotonFieldSampling::prob_eps(double eps, bool onProton, double Ein, double zIn) const {
	const double m = mass(onProton);  
	double gamma = Ein / m;
	double beta = std::sqrt(1. - 1. / gamma / gamma);
	double photonDensity = photonField->getPhotonDensity(eps, zIn);
	
	if (photonDensity != 0.) {
		double sMin = 1.1646;  // [GeV], head-on collision
		double sMax = std::max(sMin, m * m + 2. * eps / 1.e9 * Ein * (1. + beta));
		double sintegr = gaussInt([this, onProton](double s) { return this->functs(s, onProton); }, sMin, sMax);
		return photonDensity / eps / eps * sintegr / 8. / beta / Ein / Ein * 1.e18 * 1.e6;
	}
	return 0;
}


double PhotonFieldSampling::crossection(double x, bool onProton) const {
	const double m = mass(onProton); 
	const double sth = 1.1646;  // GeV^2
	const double s = m * m + 2. * m * x;
	if (s < sth)
		return 0.;
	double cross_res = 0.;
	double cross_dir = 0.;
	double cross_dir1 = 0.;
	double cross_dir2 = 0.;
	double sig_res[9];

	// first half of array: 9x proton resonance data | second half of array 9x neutron resonance data
	static const double AMRES[18] = {1.231, 1.440, 1.515, 1.525, 1.675, 1.680, 1.690, 1.895, 1.950, 1.231, 1.440, 1.515, 1.525, 1.675, 1.675, 1.690, 1.895, 1.950};
	static const double BGAMMA[18] = {5.6, 0.5, 4.6, 2.5, 1.0, 2.1, 2.0, 0.2, 1.0, 6.1, 0.3, 4.0, 2.5, 0.0, 0.2, 2.0, 0.2, 1.0};
	static const double WIDTH[18] = {0.11, 0.35, 0.11, 0.1, 0.16, 0.125, 0.29, 0.35, 0.3, 0.11, 0.35, 0.11, 0.1, 0.16, 0.150, 0.29, 0.35, 0.3};
	static const double RATIOJ[18] = {1., 0.5, 1., 0.5, 0.5, 1.5, 1., 1.5, 2., 1., 0.5, 1., 0.5, 0.5, 1.5, 1., 1.5, 2.};
	static const double AM2[2] = {0.882792, 0.880351};

	const int idx = onProton? 0 : 9;
	double SIG0[9];
	for (int i = 0; i < 9; ++i) {
		SIG0[i] = 4.893089117 / AM2[int(onProton)] * RATIOJ[i + idx] * BGAMMA[i + idx];
	}
	if (x <= 10.) {
		cross_res = breitwigner(SIG0[0], WIDTH[0 + idx], AMRES[0 + idx], x, onProton) * Ef(x, 0.152, 0.17);
		sig_res[0] = cross_res;
		for (int i = 1; i < 9; ++i) {
			sig_res[i] = breitwigner(SIG0[i], WIDTH[i + idx], AMRES[i + idx], x, onProton) * Ef(x, 0.15, 0.38);
			cross_res += sig_res[i];
		}
		// direct channel
		if ((x > 0.1) && (x < 0.6)) {
			cross_dir1 = 92.7 * Pl(x, 0.152, 0.25, 2.0)  // single pion production
					   + 40. * std::exp(-(x - 0.29) * (x - 0.29) / 0.002)
					   - 15. * std::exp(-(x - 0.37) * (x - 0.37) / 0.002);
		} else {
			cross_dir1 = 92.7 * Pl(x, 0.152, 0.25, 2.0);  // single pion production
		}
		cross_dir2 = 37.7 * Pl(x, 0.4, 0.6, 2.0);  // double pion production
		cross_dir = cross_dir1 + cross_dir2;
	}
	// fragmentation 2:
	double cross_frag2 = onProton? 80.3 : 60.2;
	cross_frag2 *= Ef(x, 0.5, 0.1) * std::pow(s, -0.34);
	// multipion production/fragmentation 1 cross section
	double cs_multidiff = 0.;
	double cs_multi = 0.;
	double cross_diffr1 = 0.;
	double cross_diffr2 = 0.;
	double cross_diffr = 0.;
	if (x > 0.85) {
		double ss1 = (x - 0.85) / 0.69;
		double ss2 = onProton? 29.3 : 26.4;
		ss2 *= std::pow(s, -0.34) + 59.3 * std::pow(s, 0.095);
		cs_multidiff = (1. - std::exp(-ss1)) * ss2;
		cs_multi = 0.89 * cs_multidiff;
		// diffractive scattering:
		cross_diffr1 = 0.099 * cs_multidiff;
		cross_diffr2 = 0.011 * cs_multidiff;
		cross_diffr = 0.11 * cs_multidiff;
		// **************************************
		ss1 = std::pow(x - 0.85, 0.75) / 0.64;
		ss2 = 74.1 * std::pow(x, -0.44) + 62. * std::pow(s, 0.08);
		double cs_tmp = 0.96 * (1. - std::exp(-ss1)) * ss2;
		cross_diffr1 = 0.14 * cs_tmp;
		cross_diffr2 = 0.013 * cs_tmp;
		double cs_delta = cross_frag2 - (cross_diffr1 + cross_diffr2 - cross_diffr);
		if (cs_delta < 0.) {
			cross_frag2 = 0.;
			cs_multi += cs_delta;
		} else {
			cross_frag2 = cs_delta;
		}
		cross_diffr = cross_diffr1 + cross_diffr2;
		cs_multidiff = cs_multi + cross_diffr;
	// in the original SOPHIA code, here is a switch for the return argument.
	// Here, only one case (compare in SOPHIA: NDIR=3) is needed.
	}
	return cross_res + cross_dir + cs_multidiff + cross_frag2;
}

double PhotonFieldSampling::Pl(double x, double xth, double xmax, double alpha) const {
	if (xth > x)
		return 0.;
	const double a = alpha * xmax / xth;
	const double prod1 = std::pow((x - xth) / (xmax - xth), a - alpha);
	const double prod2 = std::pow(x / xmax, -a);
	return prod1 * prod2;
}

double PhotonFieldSampling::Ef(double x, double th, double w) const {
	const double wth = w + th;
	if (x <= th) {
		return 0.;
	} else if ((x > th) && (x < wth)) {
		return (x - th) / w;
	} else if (x >= wth) {
		return 1.;
	} else {
		throw std::runtime_error("error in function Ef");
	}
}

double PhotonFieldSampling::breitwigner(double sigma_0, double Gamma, double DMM, double epsPrime, bool onProton) const {
	const double m = mass(onProton);
	const double s = m * m + 2. * m * epsPrime;
	const double gam2s = Gamma * Gamma * s;
	return sigma_0 * (s / epsPrime / epsPrime) * gam2s / ((s - DMM * DMM) * (s - DMM * DMM) + gam2s);
}

double PhotonFieldSampling::functs(double s, bool onProton) const {
	const double m = mass(onProton);
	const double factor = s - m * m;
	const double epsPrime = factor / 2. / m;
	const double sigmaPg = crossection(epsPrime, onProton);
	return factor * sigmaPg;
}

double PhotonFieldSampling::mass(bool onProton) const {
	const double m = onProton? 0.93827 : 0.93947;  // Gev/c^2
	return m;
}

} // namespace crpropa

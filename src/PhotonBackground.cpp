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
			this->photonEnergies.push_back(std::stod(line));
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
	return 0.1*eV;
}

void BlackbodyPhotonField::setQuantile(double q){
	if(not ((q==0.0001) or (q==0.001) or (q== 0.01)))
		throw std::runtime_error("Quantile not understood. Please use 0.01 (1%), 0.001 (0.1%) or 0.0001 (0.01%) \n");
	this->quantile = q;
}

PhotonFieldSampling::PhotonFieldSampling() {
	photonField = new CMB();
}

PhotonFieldSampling::PhotonFieldSampling(ref_ptr<PhotonField> field) {
	photonField = field;
}

double PhotonFieldSampling::sampleEps(bool onProton, double E, double z) const {
	// sample eps between epsMin ... epsMax
	double Ein = E/GeV;
	double epsMin = std::max(photonField->getMinimumPhotonEnergy(z)/eV, epsMinInteraction(onProton, Ein));
	double epsMax = photonField->getMaximumPhotonEnergy(z)/eV;
	double pEpsMax = probEpsMax(onProton, Ein, z, epsMin, epsMax);

	Random &random = Random::instance();
	for (int i = 0; i < 1000000; i++) {
		double eps = epsMin + random.rand() * (epsMax - epsMin);
		double pEps = probEps(eps, onProton, Ein, z);
		if (random.rand() * pEpsMax < pEps)
			return eps * eV;
	}
	throw std::runtime_error("error: no photon found in sampleEps, please make sure that photon field provides photons for the interaction by adapting the energy range of the tabulated photon field.");
}

double PhotonFieldSampling::epsMinInteraction(bool onProton, double Ein) const {
	// labframe energy of least energetic photon where PPP can occur
	// this kind-of ties samplingEps to the PPP and SOPHIA
	const double m = mass(onProton);
	const double Pin = sqrt(Ein * Ein - m * m);  // GeV/c
	double epsMin = 1.e9 * (1.1646 - m * m) / 2. / (Ein + Pin); // eV
	return epsMin;
}

double PhotonFieldSampling::probEpsMax(bool onProton, double Ein, double z, double epsMin, double epsMax) const {
	// find pEpsMax, which is the photon energy (eps) that we have the
	// highest (max) probability (p) to interact with
	int const nrSteps = 100;
	double pEpsMaxTested = 0.;
	double step = 0.;
	if (sampleLog){
		// sample in logspace with stepsize that is at max Δlog(E/eV) = 0.01 or otherwise dep. on size of energy range with nrSteps+1 steps log. equidis. spaced
		step = std::min(0.01, std::log10(epsMax/epsMin)/nrSteps/1.);
	} else
		step = (epsMax - epsMin) / nrSteps;

	double epsDummy = 0.;
	int i = 0;
	while (epsDummy < epsMax) {
		if (sampleLog)
			epsDummy = epsMin * pow(10,step*i);
		else
			epsDummy = epsMin + step*i;
		double p = probEps(epsDummy, onProton, Ein, z);
		if(p > pEpsMaxTested)
			pEpsMaxTested = p;
		i++;
	}
	// the following factor corrects for only trying to find the maximum on nrIteration photon energies
	// the factor should be determined in convergence tests
	double pEpsMax = pEpsMaxTested * correctionFactor;
	return pEpsMax;
}

double PhotonFieldSampling::probEps(double eps, bool onProton, double Ein, double z) const {
	// probEps returns "probability to encounter a photon of energy eps", given a primary nucleon
	// note, probEps does not return a normalized probability [0,...,1]

	double photonDensity = photonField->getPhotonDensity(eps * eV, z) * ccm / eps;
	if (photonDensity != 0.) {
		const double m = mass(onProton);
		double gamma = Ein / m;
		double beta = std::sqrt(1. - 1. / gamma / gamma);
		double sMin = 1.1646;  // [GeV^2], head-on collision
		double sMax = std::max(sMin, m * m + 2. * eps / 1.e9 * Ein * (1. + beta));
		double sintegr = gaussInt([this, onProton](double s) { return this->functs(s, onProton); }, sMin, sMax);
		return photonDensity / eps / eps * sintegr / 8. / beta / Ein / Ein * 1.e18 * 1.e6;
	}
	return 0;
}

double PhotonFieldSampling::crossection(double eps, bool onProton) const {
	const double m = mass(onProton);
	const double sth = 1.1646;  // GeV^2
	const double s = m * m + 2. * m * eps;
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
	if (eps <= 10.) {
		cross_res = breitwigner(SIG0[0], WIDTH[0 + idx], AMRES[0 + idx], eps, onProton) * Ef(eps, 0.152, 0.17);
		sig_res[0] = cross_res;
		for (int i = 1; i < 9; ++i) {
			sig_res[i] = breitwigner(SIG0[i], WIDTH[i + idx], AMRES[i + idx], eps, onProton) * Ef(eps, 0.15, 0.38);
			cross_res += sig_res[i];
		}
		// direct channel
		if ((eps > 0.1) && (eps < 0.6)) {
			cross_dir1 = 92.7 * Pl(eps, 0.152, 0.25, 2.0)  // single pion production
					   + 40. * std::exp(-(eps - 0.29) * (eps - 0.29) / 0.002)
					   - 15. * std::exp(-(eps - 0.37) * (eps - 0.37) / 0.002);
		} else {
			cross_dir1 = 92.7 * Pl(eps, 0.152, 0.25, 2.0);  // single pion production
		}
		cross_dir2 = 37.7 * Pl(eps, 0.4, 0.6, 2.0);  // double pion production
		cross_dir = cross_dir1 + cross_dir2;
	}
	// fragmentation 2:
	double cross_frag2 = onProton? 80.3 : 60.2;
	cross_frag2 *= Ef(eps, 0.5, 0.1) * std::pow(s, -0.34);
	// multipion production/fragmentation 1 cross section
	double cs_multidiff = 0.;
	double cs_multi = 0.;
	double cross_diffr1 = 0.;
	double cross_diffr2 = 0.;
	double cross_diffr = 0.;
	if (eps > 0.85) {
		double ss1 = (eps - 0.85) / 0.69;
		double ss2 = onProton? 29.3 : 26.4;
		ss2 *= std::pow(s, -0.34) + 59.3 * std::pow(s, 0.095);
		cs_multidiff = (1. - std::exp(-ss1)) * ss2;
		cs_multi = 0.89 * cs_multidiff;
		// diffractive scattering:
		cross_diffr1 = 0.099 * cs_multidiff;
		cross_diffr2 = 0.011 * cs_multidiff;
		cross_diffr = 0.11 * cs_multidiff;
		// **************************************
		ss1 = std::pow(eps - 0.85, 0.75) / 0.64;
		ss2 = 74.1 * std::pow(eps, -0.44) + 62. * std::pow(s, 0.08);
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

double PhotonFieldSampling::Pl(double eps, double epsTh, double epsMax, double alpha) const {
	if (epsTh > eps)
		return 0.;
	const double a = alpha * epsMax / epsTh;
	const double prod1 = std::pow((eps - epsTh) / (epsMax - epsTh), a - alpha);
	const double prod2 = std::pow(eps / epsMax, -a);
	return prod1 * prod2;
}

double PhotonFieldSampling::Ef(double eps, double epsTh, double w) const {
	const double wTh = w + epsTh;
	if (eps <= epsTh) {
		return 0.;
	} else if ((eps > epsTh) && (eps < wTh)) {
		return (eps - epsTh) / w;
	} else if (eps >= wTh) {
		return 1.;
	} else {
		throw std::runtime_error("error in function Ef");
	}
}

double PhotonFieldSampling::breitwigner(double sigma0, double gamma, double DMM, double epsPrime, bool onProton) const {
	const double m = mass(onProton);
	const double s = m * m + 2. * m * epsPrime;
	const double gam2s = gamma * gamma * s;
	return sigma0 * (s / epsPrime / epsPrime) * gam2s / ((s - DMM * DMM) * (s - DMM * DMM) + gam2s);
}

double PhotonFieldSampling::functs(double s, bool onProton) const {
	const double m = mass(onProton);
	const double factor = s - m * m;
	const double epsPrime = factor / 2. / m;
	const double sigmaPg = crossection(epsPrime, onProton);
	return factor * sigmaPg;
}

double PhotonFieldSampling::mass(bool onProton) const {
	const double m =  onProton ? mass_proton : mass_neutron;
	return m / GeV * c_squared;
}

void PhotonFieldSampling::setSampleLog(bool b) {
	sampleLog = b;
}

void PhotonFieldSampling::setCorrectionFactor(double factor){
	correctionFactor = factor;
}

} // namespace crpropa

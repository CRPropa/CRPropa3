#include "crpropa/PhotonBackground.h"
#include "crpropa/Common.h"
#include "crpropa/Random.h"
#include "crpropa/Units.h"

#include <vector>
#include <fstream>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <sstream>

namespace crpropa {

// Class to handle global evolution of IRB models (cf. CRPropa3-data/calc_scaling.py)
struct PhotonFieldScaling {
	bool initialized;
	std::string name;
	std::vector<double> tab_z;
	std::vector<double> tab_s;

	PhotonFieldScaling(std::string filename) {
		name = filename;
		initialized = false;
	}

	void init() {
		std::string path = getDataPath("Scaling/scaling_" + name + ".txt");
		std::ifstream infile(path.c_str());

		if (!infile.good())
			throw std::runtime_error(
					"crpropa: could not open file scaling_" + name);

		double z, s;
		while (infile.good()) {
			if (infile.peek() != '#') {
				infile >> z >> s;
				tab_z.push_back(z);
				tab_s.push_back(s);
			}
			infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}

		infile.close();
		initialized = true;
	}

	double scalingFactor(double z) {
		if (!initialized) 
#pragma omp critical(init)
			init();
		if (z > tab_z.back())
			return 0;  // zero photon background beyond maximum tabulated value
		return interpolate(z, tab_z, tab_s);
	}
};

static PhotonFieldScaling scalingKneiske04("IRB_Kneiske04");
static PhotonFieldScaling scalingStecker05("IRB_Stecker05");
static PhotonFieldScaling scalingFranceschini08("IRB_Franceschini08");
static PhotonFieldScaling scalingFinke10("IRB_Finke10");
static PhotonFieldScaling scalingDominguez11("IRB_Dominguez11");
static PhotonFieldScaling scalingGilmore12("IRB_Gilmore12");
static PhotonFieldScaling scalingStecker16_upper("IRB_Stecker16_upper");
static PhotonFieldScaling scalingStecker16_lower("IRB_Stecker16_lower");

double photonFieldScaling(PhotonField photonField, double z) {
	switch (photonField) {
	case CMB: // constant comoving photon number density
	case PF1: case PF2: case PF3: case PF4:
	case PF5: case PF6: case PF7: case PF8:
		return 1;
	case IRB:
	case IRB_Kneiske04:
		return scalingKneiske04.scalingFactor(z);
	case IRB_Stecker05:
		return scalingStecker05.scalingFactor(z);
	case IRB_Franceschini08:
		return scalingFranceschini08.scalingFactor(z);
	case IRB_Finke10:
		return scalingFinke10.scalingFactor(z);
	case IRB_Dominguez11:
		return scalingDominguez11.scalingFactor(z);
	case IRB_Gilmore12:
		return scalingGilmore12.scalingFactor(z);
	case IRB_Stecker16_upper:
		return scalingStecker16_upper.scalingFactor(z);
	case IRB_Stecker16_lower:
		return scalingStecker16_lower.scalingFactor(z);
	case URB_Protheroe96:
		if (z < 0.8) { return 1; }
		if (z < 6) { return pow((1 + 0.8) / (1 + z), 4); }
		else { return 0; }
	default:
		throw std::runtime_error("PhotonField: unknown photon background");
	}
}

std::string photonFieldName(PhotonField photonField) {
	switch (photonField) {
		case CMB: return "CMB";
		case PF1: return "PF1";
		case PF2: return "PF2";
		case PF3: return "PF3";
		case PF4: return "PF4";
		case PF5: return "PF5";
		case PF6: return "PF6";
		case PF7: return "PF7";
		case PF8: return "PF8";
		case IRB:
		case IRB_Kneiske04: return "IRB_Kneiske04";
		case IRB_Stecker05: return "IRB_Stecker05";
		case IRB_Franceschini08: return "IRB_Franceschini08";
		case IRB_Finke10: return "IRB_Finke10";
		case IRB_Dominguez11: return "IRB_Dominguez11";
		case IRB_Gilmore12: return "IRB_Gilmore12";
		case IRB_Stecker16_upper: return "IRB_Stecker16_upper";
		case IRB_Stecker16_lower: return "IRB_Stecker16_lower";
		case URB_Protheroe96: return "URB_Protheroe96";
		default:
			throw std::runtime_error("PhotonField: unknown photon background");
	}
}

CustomPhotonField::CustomPhotonField(std::string fieldPath) {
	init(fieldPath);
}

CustomPhotonField::CustomPhotonField() {
	// empty constructor for initialization in some modules
}

void CustomPhotonField::init(std::string filename) {
	std::vector< std::vector<double> > dndeps;
	std::ifstream infile(filename.c_str());
	if (!infile.good())
		throw std::runtime_error("PhotoPionProduction @ CustomPhotonField::init : could not open file " + filename);
	std::string line;
	int i = 0;
	while (std::getline(infile, line)) {
		if (line.at(0) == '#')
			continue;
		std::istringstream ss(line);
		std::vector<double> vec;
		double n;
		while (ss >> n)
			vec.push_back(n);
		if (i == 0) {
			photonEnergy = vec;
			i++;
			continue;
		}
		if (i == 1) {
			photonRedshift = vec;
			i++;
			continue;
		}
		dndeps.push_back(vec);
	}
	for (int i = 0; i < photonEnergy.size(); ++i) {
		for (int j = 0; j < photonRedshift.size(); ++j) {
			photonDensity.push_back(dndeps[i][j]);
		}
	}
	infile.close();
}

double CustomPhotonField::sampleEps(bool onProton, double Ein, double zIn) const {
/*
	- input: particle type with energy Ein at redshift z
	- output: photon energy [eV] of random photon of photon field
	- samples distribution of n(epsilon)/epsilon^2
*/
	const double zMax = photonRedshift[photonRedshift.size() - 1];
	if (zIn > zMax)
		return 0.;

	Ein /= GeV;  // SOPHIA standard unit
	// calculate pMax and its norm factor via maximum such that peps <= 1
	double cnorm = 0.;
	double pMax = 0.;
	for (int i = 0; i < photonEnergy.size(); ++i) {
		double prob = SOPHIA_probEps(photonEnergy[i], onProton, Ein, zIn);
		cnorm += prob;
		if (prob > pMax)
			pMax = prob;
	}
	pMax /= cnorm;

	// sample eps between epsMin ... epsMax
	const double epsMin = photonEnergy[0];
	const double epsMax = photonEnergy[photonEnergy.size() - 1];
	double eps;
	double peps;
	Random &random = Random::instance();
	do {
		eps = epsMin + random.rand() * (epsMax - epsMin);
		peps = SOPHIA_probEps(eps, onProton, Ein, zIn) / cnorm;
	} while (random.rand() * pMax > peps);
	return eps;
}

double CustomPhotonField::SOPHIA_probEps(double eps, bool onProton, double Ein, double zIn) const {
/*
	- input: eps [eV]
	- output: probability to encounter photon of energy eps
	- called by: sampleEps, gaussInt
*/
	const double mass = onProton? 0.93827 : 0.93947;  // Gev/c^2
	double gamma = Ein / mass;
	double beta = std::sqrt(1. - 1. / gamma / gamma);
	double photonDensity = getPhotonDensity(eps * eV, zIn) / 6.2415091e24;  // 1/(Jm³) -> 1/(eVcm³)
	if (photonDensity == 0.) {
		return 0.;
	} else {
		double sMin = 1.1646;  // head-on collision
		double sMax = std::max(sMin, mass * mass + 2. * eps / 1.e9 * Ein * (1. + beta));
		static const double x[8] = {.0950125098, .2816035507, .4580167776, .6178762444,
									.7554044083, .8656312023, .9445750230, .9894009349};
		static const double w[8] = {.1894506104, .1826034150, .1691565193, .1495959888,
									.1246289712, .0951585116, .0622535239, .0271524594};
		double xm = 0.5 * (sMax + sMin);
		double xr = 0.5 * (sMax - sMin);
		double ss = 0.;
		for (int i = 0; i < 8; ++i) {
			double dx = xr * x[i];
			ss += w[i] * (SOPHIA_functs(xm + dx, onProton) + SOPHIA_functs(xm - dx, onProton));
		}
		double sIntegral = xr * ss;
		return photonDensity / eps / eps * sIntegral / 8. / beta / Ein / Ein * 1.e24;
	}
}

double CustomPhotonField::getPhotonDensity(double eps, double z) const {
/*
	- input: photon energy, redshift
	- output: dndeps(e,z) [#/(eV cm^3)] from input file
	- called by: sampleEps
*/
	return interpolate2d(eps / eV, z, photonEnergy, photonRedshift, photonDensity) * 6.2415091e24;  // 1/(eVcm³)->1/(Jm³)
}

double CustomPhotonField::SOPHIA_crossection(double x, bool onProton) const {
/* 
	- input: photon energy [eV], specifier: 0=neutron, 1=proton
	- output: SOPHIA_crossection of nucleon-photon-interaction [area]
	- called by: SOPHIA_functs
*/  
	const double mass = onProton? 0.93827 : 0.93947;  // Gev/c^2
	double sth = 1.1646;
	double s = mass * mass + 2. * mass * x;
	if (s < sth)
		return 0.;
	// upper: proton resonance masses [GeV]
	// lower: neutron resonance masses [GeV]
	static const double AMRES[18] = {1.231, 1.440, 1.515, 1.525, 1.675, 1.680, 1.690, 1.895, 1.950,
									 1.231, 1.440, 1.515, 1.525, 1.675, 1.675, 1.690, 1.895, 1.950};
	static const double BGAMMA[18] = {5.6, 0.5, 4.6, 2.5, 1., 2.1, 2., 0.2, 1.,
									  6.1, 0.3, 4.0, 2.5, 0., 0.2, 2., 0.2, 1.};
	static const double WIDTH[18] = {0.11, 0.35, 0.11, 0.1, 0.16, 0.125, 0.29, 0.35, 0.3,
									 0.11, 0.35, 0.11, 0.1, 0.16, 0.150, 0.29, 0.35, 0.3};
	static const double RATIOJ[18] = {1., 0.5, 1., 0.5, 0.5, 1.5, 1., 1.5, 2.,
									  1., 0.5, 1., 0.5, 0.5, 1.5, 1., 1.5, 2.};
	static const double AM2[2] = {0.882792, 0.880351};

	int idx = onProton? 0 : 9;
	double SIG0[9];
	for (int i = 0; i < 9; ++i) {
		SIG0[i] = 4.893089117 / AM2[int(onProton)] * RATIOJ[i + idx] * BGAMMA[i + idx];
	}
	double cross_dir = 0.;
	double cross_res = 0.;
	if (x <= 10.) {
		cross_res = SOPHIA_breitwigner(SIG0[0], WIDTH[0 + idx], AMRES[0 + idx], x, onProton)
				  * SOPHIA_ef(x, 0.152, 0.17);
		for (int i = 1; i < 9; ++i) {
			cross_res = SOPHIA_breitwigner(SIG0[i], WIDTH[i + idx], AMRES[i + idx], x, onProton)
					   * SOPHIA_ef(x, 0.15, 0.38);
		}
		// direct channel
		double cross_dir1 = 0.;
		if ( (x > 0.1) && (x < 0.6) ) {
			cross_dir1 = 92.7 * SOPHIA_pl(x, 0.152, 0.25, 2.)  // single pion production
					   + 40.0 * std::exp(-(x - 0.29) * (x - 0.29) / 0.002)
					   - 15.0 * std::exp(-(x - 0.37) * (x - 0.37) / 0.002);
		} else {
			cross_dir1 = 92.7 * SOPHIA_pl(x, 0.152, 0.25, 2.);  // single pion production
		}
		double cross_dir2 = 37.7 * SOPHIA_pl(x, 0.4, 0.6, 2.);  // double pion production
		cross_dir = cross_dir1 + cross_dir2;
	}
	// fragmentation 2:
	double cross_frag2;
	if (onProton) {
		cross_frag2 = 80.3 * SOPHIA_ef(x, 0.5, 0.1) * std::pow(s, -0.34);
	} else {
		cross_frag2 = 60.2 * SOPHIA_ef(x, 0.5, 0.1) * std::pow(s, -0.34);
	}
	// multipion production/fragmentation 1 cross section
	double cs_multidiff = 0.;
	if (x > 0.85) {
		double ss1 = (x - 0.85) / 0.69;
		double ss2;
		if (onProton) {
			ss2 = 29.3 * std::pow(s, -0.34) + 59.3 * std::pow(s, 0.095);
		} else {
			ss2 = 26.4 * std::pow(s, -0.34) + 59.3 * std::pow(s, 0.095);
		}
		cs_multidiff = -expm1f(-ss1) * ss2;
		// diffractive scattering:
		double cross_diffr = 0.11 * cs_multidiff;
		// **************************************
		ss1 = std::pow((x - 0.85), 0.75) / 0.64;
		ss2 = 74.1 * std::pow(x, -0.44) + 62. * std::pow(s, 0.08);
		double cs_tmp = 0.96 * -expm1(-ss1) * ss2;
		double cross_diffr1 = 0.14 * cs_tmp;
		double cross_diffr2 = 0.013 * cs_tmp;
		double cs_delta = cross_frag2
						- (cross_diffr1+cross_diffr2-cross_diffr);
		double cs_multi = 0.89 * cs_multidiff;
		if (cs_delta < 0.) {
			cross_frag2 = 0.;
			cs_multi += cs_delta;
		} else {
			cross_frag2 = cs_delta;
		}
		cross_diffr = cross_diffr1 + cross_diffr2;
		cs_multidiff = cs_multi + cross_diffr;
	}
	return (cross_res + cross_dir + cs_multidiff + cross_frag2) * 1.e-34;  // µbarn to m²
}

double CustomPhotonField::SOPHIA_pl(double x, double xth, double xmax, double alpha) const {
/*
	- input: photon energy [eV], threshold [eV], max [eV], unknown [no unit]
	- output: unknown [no unit]
	- called by: SOPHIA_crossection
*/
	if (xth > x)
		return 0.;
	double a = alpha * xmax / xth;
	double prod1 = std::pow((x - xth) / (xmax - xth), (a - alpha));
	double prod2 = std::pow(x / xmax, -a);
	return prod1 * prod2;
}

double CustomPhotonField::SOPHIA_ef(double x, double th, double w) const {
/*
	- input: photon energy [eV], threshold [eV], unknown [eV]
	- output: unknown [no unit]
	- called by: SOPHIA_crossection
*/
	double wth = w + th;
	if (x <= th) {
		return 0.;
	} else if ((x > th) && (x < wth)) {
		return (x - th) / w;
	} else if (x >= wth) {
		return 1.;
	} else {
		throw std::runtime_error("error in function PhotonBackground::SOPHIA_ef");
	}
}

double CustomPhotonField::SOPHIA_breitwigner(double sigma_0, double Gamma, double DMM,
	double epsPrime, bool onProton) const {
/*
	- input: cross section [µbarn], width [GeV], mass [GeV/c^2]
	- output: Breit-Wigner SOPHIA_crossection of a resonance widh width Gamma
	- called by: SOPHIA_crossection
*/
	const double mass = onProton? 0.93827 : 0.93947;  // Gev/c^2
	double s = mass * mass + 2. * mass * epsPrime;
	double gam2s = Gamma * Gamma * s;
	return sigma_0 * (s / epsPrime / epsPrime) * gam2s
				   / ((s - DMM*DMM) * (s - DMM * DMM) + gam2s);
}

double CustomPhotonField::SOPHIA_functs(double s, bool onProton) const {
/*
	- input: s [GeV^2] 
	- output: (s-p^2)*sigma_(nucleon/gamma) [GeV^2*area]
	- called by: sampleEps
*/  
	const double mass = onProton? 0.93827 : 0.93947;  // Gev/c^2
	double factor = s - mass * mass;
	double epsPrime = factor / 2. / mass;
	double sigma_pg = SOPHIA_crossection(epsPrime, onProton) / 1.e-34;  // m² to µbarn
	return factor * sigma_pg;
}

} // namespace crpropa

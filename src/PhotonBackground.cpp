#include "crpropa/PhotonBackground.h"
#include "crpropa/Common.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <vector>
#include <fstream>
#include <stdexcept>
#include <limits>
#include <cmath>

namespace crpropa {

// Class to handle global evolution of IRB models (cf. CRPropa3-data/calc_scaling.py)
struct PhotonFieldScaling {
	bool isInitialized;
	std::string name;
	std::vector<double> tab_z;
	std::vector<double> tab_s;

	PhotonFieldScaling(std::string filename) {
		name = filename;
		isInitialized = false;
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
		isInitialized = true;
	}

	double scalingFactor(double z) {
		if (!isInitialized) 
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
	case CMB:
		return 1;  // constant comoving photon number density
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
		if (z < 0.8)
			return 1;
		if (z < 6)
			return pow((1 + 0.8) / (1 + z), 4);
		else
			return 0;
	default:
		throw std::runtime_error("PhotonField: unknown photon background");
	}
}

std::string photonFieldName(PhotonField photonField) {
	switch (photonField) {
	case CMB:
		return "CMB";
	case IRB:
	case IRB_Kneiske04:
		return "IRB_Kneiske04";
	case IRB_Stecker05:
		return "IRB_Stecker05";
	case IRB_Franceschini08:
		return "IRB_Franceschini08";
	case IRB_Finke10:
		return "IRB_Finke10";
	case IRB_Dominguez11:
		return "IRB_Dominguez11";
	case IRB_Gilmore12:
		return "IRB_Gilmore12";
	case IRB_Stecker16_upper:
		return "IRB_Stecker16_upper";
	case IRB_Stecker16_lower:
		return "IRB_Stecker16_lower";
	case URB_Protheroe96:
		return "URB_Protheroe96";
	default:
		throw std::runtime_error("PhotonField: unknown photon background");
	}
}

PhotonFieldSampling::PhotonFieldSampling() {
	bgFlag = 0;
}

PhotonFieldSampling::PhotonFieldSampling(int flag) {
	if (flag != 1 && flag != 2)
		throw std::runtime_error("error: incorrect background flag. Must be 1 (CMB) or 2 (IRB_Kneiske04).");
	bgFlag = flag;
}

double PhotonFieldSampling::sample_eps(bool onProton, double E_in, double z_in) const {
	if (bgFlag == 0)
		throw std::runtime_error("error: select photon field first: 1 (CMB) or 2 (IRB_Kneiske04)");

	const double mass = onProton? 0.93827 : 0.93947;  // Gev/c^2
	const double P_in = sqrt(E_in * E_in - mass * mass);  // GeV/c

	double eps = 0.;
	if (bgFlag == 1) {
		// CMB
		const double tbb = 2.73 * (1. + z_in);
		const double epsMin = (1.1646 - mass * mass) / 2. / (E_in + P_in) * 1.e9;
		const double epsMax = 0.007 * tbb;
		if (epsMin > epsMax) {
			std::cout << "sample_eps (CMB): CMF energy is below threshold for nucleon energy " << E_in << " GeV !" << std::endl;
			return 0.;
		}

		const double cnorm = gaussInt("prob_eps", epsMin, epsMax, onProton, E_in, z_in);
		const double epskt = 8.619e-5 * tbb;
		const double epspmax = (3.e-3 * std::pow(E_in * epskt * 1.e-9, -0.97) + 0.047) / 3.9e2 * tbb;
		const double pmaxc = prob_eps(epspmax, onProton, E_in, z_in) / cnorm;
		const double facpmax = 1.6;
		const double pMax = facpmax * pmaxc;

		// sample eps between epsMin ... epsMax
		double peps = 0.;
		Random &random = Random::instance();
		do {
			eps = epsMin + random.rand() * (epsMax - epsMin);
			peps = prob_eps(eps, onProton, E_in, z_in) / cnorm;
		} while (random.rand() * pMax > peps);
	}

	if (bgFlag == 2) {
		// IRB_Kneiske04     
		const double epsMin = std::max(0.00395, 1.e9 * (1.1646 - mass * mass) / 2. / (E_in + P_in));  // eV
		const double epsMax = 12.2;  // eV
		if (epsMin > epsMax) {
			std::cout << "sample_eps (IRB): CMF energy is below threshold for nucleon energy " << E_in << " GeV !" << std::endl;
			return 0.;
		}
		const int i_max = static_cast<int>(10. * std::log(epsMax / epsMin)) + 1;
		const double de = std::log(epsMax / epsMin) / i_max;
		double rmax = 0.;
		double eps_dum = 0.;
		double dum = 0.;
		for (int i = 0; i < i_max; ++i) {
			eps_dum = epsMin * std::exp(i * de);
			dum = eps_dum * eps_dum * getPhotonDensity(eps_dum, z_in);
			if (dum > rmax)
				rmax = dum;
		}
		const double beta = 4.;
		const double e1 = std::pow(epsMin, 1. - beta);
		const double e2 = std::pow(epsMax, 1. - beta);
		bool keepTrying = true;
		int i_rep = 0;
		Random &random = Random::instance();
		do {
			if (i_rep >= 100000) {
				keepTrying = false;
				return 0.;
			}
			eps = std::pow(random.rand() * (e1 - e2) + e2, 1./(1. - beta));
			i_rep++;
			if (random.rand() < eps * eps * getPhotonDensity(eps, z_in) / rmax)
				keepTrying = false;
		} while (keepTrying);
	}
	return eps * eV;
}

double PhotonFieldSampling::prob_eps(double eps, bool onProton, double E_in, double z_in) const {
	const double mass = onProton? 0.93827 : 0.93947;  // Gev/c^2
	double gamma = E_in / mass;
	double beta = std::sqrt(1. - 1. / gamma / gamma);
	double photonDensity = getPhotonDensity(eps, z_in);
	if (photonDensity == 0.) {
		return 0.;
	} else {
		double smin = 1.1646;  // [GeV], head-on collision
		double smax = std::max(smin, mass * mass + 2. * eps / 1.e9 * E_in * (1. + beta));
		double sintegr = gaussInt("functs", smin, smax, onProton, E_in, z_in);
		return photonDensity / eps / eps * sintegr / 8. / beta / E_in / E_in * 1.e18 * 1.e6;
	}
}

double PhotonFieldSampling::getPhotonDensity(double eps, double z_in) const {
	if (bgFlag == 1) {
		// CMB
		return 1.318e13 * eps * eps / (std::exp(eps / (8.619e-5 * 2.73)) - 1.);
	}

	if (bgFlag == 2) {
		// IR background from Primack et al. (1999) 
		const double ZMAX_IR = 5.;
		if (z_in > ZMAX_IR)
			return 0.;
		const double X = 1.2398 * (1. + z_in) / eps;
		if (X > 500.)
			return 0.;

		static const double XData[15] = {-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5};  // log10(eV)
		static const double YData[15] = {-0.214401, 0.349313, 0.720354, 0.890389, 1.16042, 1.24692, 1.06525, 0.668659, 0.536312, 0.595859, 0.457456, 0.623521, 1.20208, 1.33657, 1.04461};  // log10(nW/m^2/sr)
		
		// interpolate
		if (std::log10(X) <= XData[0])
			return 0.;
		if (std::log10(X) >= XData[14]) {
			double result = (YData[14] - YData[13]) / (XData[14] - XData[13]) * (std::log10(X) - XData[13]) + YData[13];
			return std::pow(10., result);
		}
		int index = 1;
		do {
			index++;
		} while (XData[index] < std::log10(X));
		double result = (YData[index] - YData[index - 1]) / (XData[index] - XData[index - 1]) * (std::log10(X) - XData[index - 1]) + YData[index - 1];
		result = std::pow(10., result);

		const double fluxConversion = 3.82182e3;  // conversion from nW/cm^3/sr to eV/cm^3
		return result * std::pow(1. + z_in, 4.) / (eps * eps) / fluxConversion;
	}
}

double PhotonFieldSampling::crossection(double x, bool onProton) const {
	const double mass = onProton? 0.93827 : 0.93947;  // Gev/c^2
	const double sth = 1.1646;  // GeV^2
	const double s = mass * mass + 2. * mass * x;
	if (s < sth)
		return 0.;
	double cross_res = 0.;
	double cross_dir = 0.;
	double cross_dir1 = 0.;
	double cross_dir2 = 0.;
	double sig_res[9];

	// upper half of array: 9x proton resonance data | lower half of array 9x neutron resonance data
	static const double AMRES[18] = {1.231, 1.440, 1.515, 1.525, 1.675, 1.680, 1.690, 1.895, 1.950,
							  1.231, 1.440, 1.515, 1.525, 1.675, 1.675, 1.690, 1.895, 1.950};
	static const double BGAMMA[18] = {5.6, 0.5, 4.6, 2.5, 1.0, 2.1, 2.0, 0.2, 1.0,
							   6.1, 0.3, 4.0, 2.5, 0.0, 0.2, 2.0, 0.2, 1.0};
	static const double WIDTH[18] = {0.11, 0.35, 0.11, 0.1, 0.16, 0.125, 0.29, 0.35, 0.3,
							  0.11, 0.35, 0.11, 0.1, 0.16, 0.150, 0.29, 0.35, 0.3};
	static const double RATIOJ[18] = {1., 0.5, 1., 0.5, 0.5, 1.5, 1., 1.5, 2.,
							   1., 0.5, 1., 0.5, 0.5, 1.5, 1., 1.5, 2.};
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
	double cross_frag2 = 0.;
	if (onProton) {
		cross_frag2 = 80.3 * Ef(x, 0.5, 0.1) * std::pow(s, -0.34);
	} else {
		cross_frag2 = 60.2 * Ef(x, 0.5, 0.1) * std::pow(s, -0.34);
	}
	// multipion production/fragmentation 1 cross section
	double cs_multidiff = 0.;
	double cs_multi = 0.;
	double cross_diffr1 = 0.;
	double cross_diffr2 = 0.;
	double cross_diffr = 0.;
	if (x > 0.85) {
		double ss1 = (x - 0.85) / 0.69;
		double ss2 = 0.;
		if (onProton) {
			ss2 = 29.3 * std::pow(s, -0.34) + 59.3 * std::pow(s, 0.095);
		} else {
			ss2 = 26.4 * std::pow(s, -0.34) + 59.3 * std::pow(s, 0.095);
		}
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
		double cs_delta = cross_frag2
						- (cross_diffr1+cross_diffr2-cross_diffr);
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
	const double mass = onProton? 0.93827 : 0.93947;  // Gev/c^2
	const double s = mass * mass + 2. * mass * epsPrime;
	const double gam2s = Gamma * Gamma * s;
	return sigma_0 * (s / epsPrime / epsPrime) * gam2s / ((s - DMM * DMM) * (s - DMM * DMM) + gam2s);
}


double PhotonFieldSampling::functs(double s, bool onProton) const {
	const double mass = onProton? 0.93827 : 0.93947;  // Gev/c^2
	const double factor = s - mass * mass;
	const double epsPrime = factor / 2. / mass;
	const double sigma_pg = crossection(epsPrime, onProton);
	return factor * sigma_pg;
}

double PhotonFieldSampling::gaussInt(std::string type, double A, double B, bool onProton, double E_in, double z_in) const {
	static const double X[8] = {.0950125098, .2816035507, .4580167776, .6178762444, .7554044083, .8656312023, .9445750230, .9894009349};
	static const double W[8] = {.1894506104, .1826034150, .1691565193, .1495959888, .1246289712, .0951585116, .0622535239, .0271524594};
	const double XM = 0.5 * (B + A);
	const double XR = 0.5 * (B - A);
	double SS = 0.;
	for (int i = 0; i < 8; ++i) {
		double DX = XR * X[i];
		if (type == "functs") {
			SS += W[i] * (functs(XM + DX, onProton) + functs(XM - DX, onProton));
		} else if (type == "prob_eps") {
			SS += W[i] * (prob_eps(XM + DX, onProton, E_in, z_in) + prob_eps(XM - DX, onProton, E_in, z_in));
		} else {
			throw std::runtime_error("gaussInt: type incorrectly specified");
		}
	}
	return XR * SS;
}

} // namespace crpropa

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


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// custom photon field methods related to SAMPLING
// These methods are taken from SOPHIA.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Photon_Field::Photon_Field(std::string fieldPath) {
    // constructor
    init(fieldPath);
}


Photon_Field::Photon_Field() {
    // constructor 2
}


double Photon_Field::sample_eps(bool onProton, double E_in, double z_in) const {
/*
    - input: particle type 0=neutron, 1=proton, its energy [GeV], its redshift
    - output: photon energy [eV] of random photon of photon field
    - samples distribution of n(epsilon)/epsilon^2
*/ 
    const double mass = onProton? 0.93827 : 0.93947;  // Gev/c^2
    const double z_max = redshift[redshift.size()-1];
    if (z_in > z_max) { return 0.; }
    double z_pos;
    double smallestDiff = z_max;
    for (int i = 0; i < redshift.size(); ++i) {
        double diff = std::abs(z_in-redshift[i]);
        if (diff < smallestDiff) {
            smallestDiff = diff;
            z_pos = i;
        }
    }

    const double epsMin = energy[0];
    const double epsMax = energy[energy.size()-1];
    double eps;
    double p = std::sqrt(E_in*E_in-mass*mass);
    double peps;
    double cnorm = gaussInt("prob_eps", epsMin, epsMax, onProton, E_in, z_pos);

    // calculate pMax
    double highest = 0.;
    int closestPos;
    for (int i = 0; i < energy.size(); ++i) {
        if (dn_deps[i][z_pos] > highest) {
            highest = dn_deps[i][z_pos];
            closestPos = i;
        }
    }
    double eps_pMax = energy[closestPos];
    double pMax = prob_eps(eps_pMax, onProton, E_in, z_pos)/cnorm;
    if ( (pMax < 0.01) || (pMax > 1.) ) { pMax = 1.; }

    // sample eps randomly between epsMin ... epsMax
    Random &random = Random::instance();
    do {
        eps = epsMin+random.rand()*(epsMax-epsMin);
        peps = prob_eps(eps, onProton, E_in, z_pos)/cnorm;
    } while (random.rand()*pMax > peps);
    return eps;
}  // end sample_eps


void Photon_Field::init(std::string filename) {
    std::ifstream infile(filename.c_str());
    if (!infile.good()) {
        throw std::runtime_error("PhotoPionProduction @ Photon_Field::init : could not open file " + filename);
    }
    std::string line;
    int i = 0;
    while ( std::getline(infile, line) ) {
        if (line.find('#') == 0 ) continue;
        std::istringstream ss(line);
        std::vector<double> vec;
        double n;
        while (ss >> n) vec.push_back(n);
        if (i == 0) { energy = vec; i++; continue; }
        if (i == 1) { redshift = vec; i++; continue; }
        dn_deps.push_back(vec);
    }
    infile.close();
}  // end init


double Photon_Field::prob_eps(double eps, bool onProton, double E_in, int z_pos) const {
/*
    - input: eps [eV]
    - output: probability to encounter photon of energy eps
    - called by: sample_eps, gaussInt
*/
    const double mass = onProton? 0.93827 : 0.93947;  // Gev/c^2
    double gamma = E_in/mass;
    double beta = std::sqrt(1.-1./gamma/gamma);
    double photonDensity = get_photonDensity(eps, z_pos);
    if (photonDensity == 0.) {
        return 0.;
    } else {
        double smin = 1.1646;  // head-on collision
        double smax = std::max(smin, mass*mass+2.*eps/1.e9*E_in*(1.+beta));
        double sintegr = gaussInt("functs", smin, smax, onProton, E_in, z_pos);
        return photonDensity/eps/eps*sintegr/8./beta/E_in/E_in*1.e18*1.e6;
    }
}  // end prob_eps


double Photon_Field::get_photonDensity(double eps, int z_pos) const {
/*
    - input: photon energy [eV], redshift
    - output: dn_deps(e,z) [#/(eV cm^3)] from input file
    - called by: sample_eps
*/
    // find closest dn_deps
    double smallestDiff = energy[energy.size()-1];
    int closestPos;
    for (int i = 0; i < energy.size(); ++i) {
        double diff = std::abs(eps-energy[i]);
        if (diff < smallestDiff) {
            smallestDiff = diff;
            closestPos = i;
        }
    }
    // linear interpolation of energy
    double realDiff = eps-energy[closestPos];
    double rho;
    if (realDiff >= 0.) {
        rho = realDiff / (energy[closestPos+1] - energy[closestPos]);
        return (1.-rho)*dn_deps[closestPos][z_pos]
               + rho*dn_deps[closestPos+1][z_pos];
    } else {
        rho = 1. - (std::abs(realDiff)/energy[closestPos-1]);
        return (1.-rho)*dn_deps[closestPos-1][z_pos]
               + rho*dn_deps[closestPos][z_pos];
    }
}  // end get_photonDensity


double Photon_Field::crossection(double x, bool onProton) const {
/* 
    - input: photon energy [eV], specifier: 0=neutron, 1=proton
    - output: crossection of nucleon-photon-interaction [area]
    - called by: functs
*/  
    const double mass = onProton? 0.93827 : 0.93947;  // Gev/c^2
    double sth = 1.1646;
    double s = mass*mass + 2.*mass*x;
    if (s < sth) { return 0.; }
    double cross_res = 0.;
    double cross_dir = 0.;
    double cross_dir1 = 0.;
    double cross_dir2 = 0.;
    double sig_res[9];

    // upper: ppppppppp
    // lower: nnnnnnnnn
    const double AMRES[18] = { 1.231, 1.440, 1.515, 1.525, 1.675, 1.680, 1.690, 1.895, 1.950,
                               1.231, 1.440, 1.515, 1.525, 1.675, 1.675, 1.690, 1.895, 1.950 };
    const double BGAMMA[18] = { 5.6, .5, 4.6, 2.5, 1., 2.1, 2., .2, 1.,
                                6.1, .3, 4., 2.5, 0.0, .2, 2., .2, 1.};
    const double WIDTH[18] = { .11, .35, .11, .1, .16, .125, .29, .35, .3,
                               .11, .35, .11, .1, .16, .15, .29, .35, .3};
    const double RATIOJ[18] = { 1., .5, 1., .5, .5, 1.5, 1., 1.5, 2.,
                                1., .5, 1., .5, .5, 1.5, 1., 1.5, 2.};
    const double AM2[2] = { 0.882792, 0.880351 };

    // int idx = (onProton == 0)? 9:0;  // neutron = 0, proton = 1
    int idx = onProton? 0:9;  // neutron = 0, proton = 1
    double SIG0[9];
    for (int i = 0; i < 9; ++i) {
        SIG0[i] = 4.893089117/AM2[int(onProton)]*RATIOJ[i+idx]*BGAMMA[i+idx];
    }

    if (x <= 10.) {
        cross_res = breitwigner(SIG0[0], WIDTH[0+idx], AMRES[0+idx], x, onProton)
                  * Ef(x, 0.152, 0.17);
        sig_res[0] = cross_res;
        for (int i = 1; i < 9; ++i) {
            sig_res[i] = breitwigner(SIG0[i], WIDTH[i+idx], AMRES[i+idx], x, onProton)
                       * Ef(x, 0.15, 0.38);
            cross_res += sig_res[i];
        }
        // direct channel
        if ( (x > 0.1) && (x < 0.6) ) {
            cross_dir1 = singleback(x)
                       + 40.*std::exp(-(x-0.29)*(x-0.29) / 0.002)
                       - 15.*std::exp(-(x-0.37)*(x-0.37) / 0.002);
        } else {
            cross_dir1 = singleback(x);
        }
        cross_dir2 = twoback(x);
        cross_dir = cross_dir1 + cross_dir2;
    }
    // fragmentation 2:
    double cross_frag2;
    if (onProton) {
        cross_frag2 = 80.3*Ef(x, 0.5, 0.1) * std::pow(s, -0.34);
    } else {
        cross_frag2 = 60.2*Ef(x, 0.5, 0.1) * std::pow(s, -0.34);
    }
    // multipion production/fragmentation 1 cross section
    double cs_multidiff;
    double cs_multi;
    double cross_diffr1;
    double cross_diffr2;
    double cross_diffr;
    if (x > 0.85) {
        double ss1 = (x-0.85)/.69;
        double ss2;
        if (onProton) {
            ss2 = 29.3*std::pow(s, -0.34) + 59.3*std::pow(s, 0.095);
        } else {
            ss2 = 26.4*std::pow(s, -0.34) + 59.3*std::pow(s, 0.095);
        }
        cs_multidiff = (1.-std::exp(-ss1))*ss2;
        cs_multi = 0.89*cs_multidiff;
        // diffractive scattering:
        cross_diffr1 = .099*cs_multidiff;
        cross_diffr2 = .011*cs_multidiff;
        cross_diffr = .11*cs_multidiff;
        // **************************************
        ss1 = std::pow((x-.85), .75)/.64;
        ss2 = 74.1*std::pow(x, -0.44) + 62.*std::pow(s, 0.08);
        double cs_tmp = 0.96*(1.-std::exp(-ss1))*ss2;
        cross_diffr1 = 0.14*cs_tmp;
        cross_diffr2 = 0.013*cs_tmp;
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
    // *****************************************
    } else {
        cross_diffr = 0.;
        cross_diffr1 = 0.;
        cross_diffr2 = 0.;
        cs_multidiff = 0.;
        cs_multi = 0.;
    }
    // in the original SOPHIA code, this is a switch.
    // Here, only one case (NDIR=3) is needed.
    return cross_res+cross_dir+cs_multidiff+cross_frag2;
}  // end crossection


double Photon_Field::Pl(double x, double xth, double xmax, double alpha) const {
/*
    - input: photon energy [eV], unknown, unknown, unknown
    - output: unknown.
    - called by: crossection
*/
    if (xth > x) { return 0.; }
    double a = alpha*xmax/xth;
    double prod1 = std::pow((x-xth)/(xmax-xth), (a-alpha));
    double prod2 = std::pow(x/xmax, -a);
    return prod1*prod2;
}  // end Pl


double Photon_Field::Ef(double x, double th, double w) const {
/*
    - input: photon energy [eV], unknown, unknown
    - output: unknown
    - called by: crossection
*/
    double wth = w+th;
    if (x <= th) {
        return 0.;
    } else if ( (x > th) && (x < wth) ) {
        return (x-th)/w;
    } else if (x >= wth) {
        return 1.;
    } else {
        throw std::runtime_error("error in function Ef");
    }
}  // end Ef


double Photon_Field::breitwigner(double sigma_0, double Gamma, double DMM, double eps_prime, bool onProton) const {
/*
    - input: cross section [µbarn], width [GeV], mass [GeV/c^2]
    - output: Breit-Wigner crossection of a resonance widh width Gamma
    - called by: crossection
*/
    const double mass = onProton? 0.93827 : 0.93947;  // Gev/c^2
    double s = mass*mass + 2.*mass*eps_prime;
    double gam2s = Gamma*Gamma*s;
    return sigma_0 * (s/eps_prime/eps_prime)*gam2s
                   / ((s-DMM*DMM)*(s-DMM*DMM)+gam2s);
}  // end breitwigner


double Photon_Field::singleback(double x) const {
/*
    - single pion channel
    - called by: crossection
*/    
    return 92.7 * Pl(x, .152, .25, 2.);
}  // end singleback


double Photon_Field::twoback(double x) const {
/*
    - two pion production
    - called by: crossection
*/
    return 37.7 * Pl(x, .4, .6, 2.);
}  // end twoback


double Photon_Field::functs(double s, bool onProton) const {
/*
    - input: s [GeV^2] 
    - output: (s-p^2)*sigma_(nucleon/gamma) [GeV^2*area]
    - called by: sample_eps
*/  
    const double mass = onProton? 0.93827 : 0.93947;  // Gev/c^2
    double factor = s - mass*mass;
    double epsprime = factor/2./mass;
    double sigma_pg = crossection(epsprime, onProton);
    return factor*sigma_pg;
}  // end functs


double Photon_Field::gaussInt(std::string type, double A, double B, bool onProton, double E_in, int z_pos) const {
/* 
    - input:  type: specifier of function over which to integrate,
              integration limits A and B
    - output: 8-points gauss-Legendre integral
    - called by: sample_eps, prob_eps
*/
    const double X[8] = { .0950125098, .2816035507, .4580167776, .6178762444,
                          .7554044083, .8656312023, .9445750230, .9894009349 };
    const double W[8] = { .1894506104, .1826034150, .1691565193, .1495959888,
                          .1246289712, .0951585116, .0622535239, .0271524594 };
    double XM = 0.5*(B+A);
    double XR = 0.5*(B-A);
    double SS = 0.;
    for (int i = 0; i < 8; ++i) {
        double DX = XR*X[i];
        if (type == "functs") {
            SS += W[i] * (functs(XM+DX, onProton) + functs(XM-DX, onProton));
        } else if (type == "prob_eps") {
            SS += W[i] * (prob_eps(XM+DX, onProton, E_in, z_pos) + prob_eps(XM-DX, onProton, E_in, z_pos));
        } else {
            throw std::runtime_error("gaussInt: type incorrectly specified");
        }
    }
    return XR*SS;
}  // end gaussInt

} // namespace crpropa

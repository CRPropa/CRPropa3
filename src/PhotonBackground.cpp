#include "crpropa/PhotonBackground.h"
#include "crpropa/Common.h"

#include <vector>
#include <fstream>
#include <stdexcept>
#include <limits>
#include <math.h>

namespace crpropa {

// Class to handle global evolution of IRB models (cf. CRPropa3-data/calc_scaling.py)
struct PhotonFieldScaling {
	std::vector<double> tab_z;
	std::vector<double> tab_s;

	PhotonFieldScaling(std::string filename) {
		std::string path = getDataPath("Scaling/scaling_" + filename + ".txt");
		std::ifstream infile(path.c_str());

		if (!infile.good())
			throw std::runtime_error(
					"crpropa: could not open file " + filename);

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
	}

	double scalingFactor(double z) {
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

} // namespace crpropa

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
		std::string path = getDataPath(filename);
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

static PhotonFieldScaling scalingKneiske04("scaling_Kneiske04.txt");
static PhotonFieldScaling scalingStecker05("scaling_Stecker05.txt");
static PhotonFieldScaling scalingFranceschini08("scaling_Franceschini08.txt");
static PhotonFieldScaling scalingFinke10("scaling_Finke10.txt");
static PhotonFieldScaling scalingDominguez11("scaling_Dominguez11.txt");
static PhotonFieldScaling scalingGilmore12("scaling_Gilmore12.txt");

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

} // namespace crpropa

#include "crpropa/PhotonBackground.h"
#include "crpropa/Common.h"

#include <vector>
#include <stdexcept>

namespace crpropa {

// Overall redshift scaling of the Kneiske et al. 2004 IRB, astro-ph/0309141
// The scaling is calculated in data-tools/PhotonField/Kneiske2004_IRB_scaling.py
double a[10] = { -1, 0, 0.2, 0.4, 0.6, 1, 2, 3, 4, 5 };
static std::vector<double> zKneiske(a, a + sizeof(a) / sizeof(double));
double b[10] = { 1.250665, 1., 0.9749867, 0.93999977, 0.88430409, 0.64952017,
		0.27170436, 0.130244, 0.05971749, 0. };
static std::vector<double> sKneiske(b, b + sizeof(b) / sizeof(double));

double photonFieldScaling(PhotonField photonField, double z) {
	switch (photonField) {
	case CMB:
		return 1; // CMB-like scaling
	case IRB:
	case IRB_Kneiske04:
	case IRB_Stecker05:
	case IRB_Franceschini08:
	case IRB_Finke10:
	case IRB_Dominguez11:
	case IRB_Gilmore12:
		return interpolate(z, zKneiske, sKneiske);
	case IRB_withRedshift_Kneiske04:
	case IRB_withRedshift_Franceschini08:
	case IRB_withRedshift_Finke10:
	case IRB_withRedshift_Gilmore12:
		return 1;
	default:
		throw std::runtime_error("PhotonField: unknown photon background");
	}
}

} // namespace crpropa

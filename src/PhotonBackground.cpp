#include "crpropa/PhotonBackground.h"
#include "crpropa/Common.h"

#include <vector>
#include <math.h>

namespace crpropa {

// Overall redshift scaling of the Kneiske IRB (see data/PhotonField/KneiskeIRB.py and Kneiske et al. 2004, astro-ph/0309141)
double a[9] = { 0, 0.2, 0.4, 0.6, 1, 2, 3, 4, 5 };
static std::vector<double> zKneiske(a, a + sizeof(a) / sizeof(double));
double b[9] = { 1, 1.6937, 2.5885, 3.6178, 5.1980, 7.3871, 8.5471, 7.8605, 0 };
static std::vector<double> sKneiske(b, b + sizeof(b) / sizeof(double));

double photonFieldScaling(int photonField, double z) {
	if (photonField == IRB)
		return interpolate(z, zKneiske, sKneiske);
	return pow(1 + z, 3); // CMB-like scaling
}

} // namespace crpropa

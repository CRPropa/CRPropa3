#include "EleCa/Common.h"
#include "EleCa/Constants.h"

#include <string>
#include <vector>
#include <cstdlib>

namespace eleca {

double z2Mpc(double z) {
	if (z < 0.4) {
		return ((double) C_speed / 1000. / H0) * z;
	} else {
		// AV Uryson, Physics Particles and Nuclei, 2006, Vol. 37, No. 3, pp. 347   67
		// Assuming flat-matter-dominated cosmology. Error is negligible
		return ((double) C_speed / 1000. / H0) * (2. / 3.)
				* (1 - pow(1. + z, -1.5));
	}
}

double Mpc2z(double D) {
	if (D < 1700.) {
		return (double) D / ((double) C_speed / 1000. / H0);
	} else {
		// AV Uryson, Physics Particles and Nuclei, 2006, Vol. 37, No. 3, pp. 347   67
		// Assuming flat-matter-dominated cosmology. Error is negligible
		return pow(1 - (double) D / ((2. / 3.) * (double) C_speed / 1000. / H0),
				-2. / 3.) - 1;
	}
}

static double (*_Uniform)(double min, double max) = 0;

void setUniformCallback(double (*uniform)(double min, double max)) {
	_Uniform = uniform;
}

double Uniform(double min, double max) {
	if (_Uniform)
		return _Uniform(min, max);
	else
		//return crpropa::Random::instance().randUniform(min, max);
		return min + (max - min) * ::drand48();
}

} // namespace eleca


#include "mpc/magneticField/sphMagneticFieldGrid.hpp"
#include "mpc/Units.h"
//#include "gadget/Vector3.hpp"
//#include "gadget/SmoothParticle.hpp"
//#include "gadget/MagneticField.hpp"

namespace mpc {

SPH_MagneticFieldGrid::SPH_MagneticFieldGrid(size_t n, double spacing,
		Vector3 origin, std::string filename) :
		MagneticFieldGrid(n, spacing, origin) {
	this->filename = filename;
	initialize();
}

void SPH_MagneticFieldGrid::initialize() {
	// here be the initialization using gadget
}

} // namespace mpc


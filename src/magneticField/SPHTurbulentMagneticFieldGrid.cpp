#include "mpc/magneticField/SPHTurbulentMagneticFieldGrid.h"

namespace mpc {

void SPHTurbulentMagneticFieldGrid::setModulation(std::string filename,
		Vector3d origin, double size, size_t samples, double exp, double norm) {
	exponent = exp;
	normalization = norm;
	sphField.reset(new SPHMagneticField(origin, size, samples, filename));
}

Vector3d SPHTurbulentMagneticFieldGrid::getField(
		const Vector3d &position) const {
	Vector3d b = MagneticFieldGrid::getField(position);
	double rho = sphField->getRho(position);
	return b * pow(rho, exponent) * normalization;
}

double SPHTurbulentMagneticFieldGrid::getRho(const Vector3d &position) const {
	return sphField->getRho(position);
}

} // namespace mpc

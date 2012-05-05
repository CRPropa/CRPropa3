#include "mpc/magneticField/SPHMagneticField.h"
#include "mpc/Units.h"

#include "kiss/logger.h"

namespace mpc {

SPHMagneticField::SPHMagneticField(Vector3d origin, double size, size_t gridSize,
		const std::string filename) :
		field(gridSize) {
	database.open(filename);
	gadget::Vector3f v = gadget::Vector3f(origin.x, origin.y, origin.z)
			/ kpc;
	field.init(v, size / kpc, database);
}

SPHMagneticField::SPHMagneticField(size_t gridSize, const std::string filename) :
		field(gridSize) {
	database.open(filename);
}

Vector3d SPHMagneticField::getField(const Vector3d &position) const {
	gadget::Vector3f r = gadget::Vector3f(position.x, position.y,
			position.z);
	gadget::Vector3f b = field.getField(r / kpc);
	Vector3d bField = Vector3d(b.x, b.y, b.z) * gauss;
	return bField;
}

void SPHMagneticField::updateSimulationVolume(const Vector3d &origin,
		double size) {
	gadget::Vector3f v = gadget::Vector3f(origin.x, origin.y, origin.z)
			/ kpc;
	field.init(v, size / kpc, database);
}

SPHMagneticFieldGrid::SPHMagneticFieldGrid(Vector3d origin, double size,
		size_t samples, const std::string filename) :
		field(samples) {
	database.open(filename);
	gadget::Vector3f v = gadget::Vector3f(origin.x, origin.y, origin.z)
			/ kpc;
	field.init(v, size / kpc, database);
}

SPHMagneticFieldGrid::SPHMagneticFieldGrid(size_t samples,
		const std::string filename) :
		field(samples) {
	database.open(filename);
}

Vector3d SPHMagneticFieldGrid::getField(const Vector3d &position) const {
	gadget::Vector3f r = gadget::Vector3f(position.x, position.y,
			position.z);
	gadget::Vector3f b = field.getField(r / kpc);
	Vector3d bField = Vector3d(b.x, b.y, b.z) * gauss;
	return bField;
}

void SPHMagneticFieldGrid::updateSimulationVolume(const Vector3d &origin,
		double size) {
	gadget::Vector3f v = gadget::Vector3f(origin.x, origin.y, origin.z)
			/ kpc;
	field.init(v, size / kpc, database);
}

} // namespace mpc


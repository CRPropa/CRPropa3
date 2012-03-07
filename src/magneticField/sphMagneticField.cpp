#include "mpc/magneticField/sphMagneticField.h"
#include "mpc/Units.h"

#include "kiss/logger.h"

namespace mpc {

SPHMagneticField::SPHMagneticField(Vector3 origin, double size, size_t gridSize,
		const std::string filename) :
		field(gridSize) {
	database.open(filename);
	gadget::Vector3f v = gadget::Vector3f(origin.x(), origin.y(), origin.z())
			/ kpc;
	field.init(v, size / kpc, database);
	KISS_LOG_INFO
		<< "mpc::SPHMagneticField: ctor" << std::endl;
}

Vector3 SPHMagneticField::getField(const Vector3 &position) const {
	gadget::Vector3f r = gadget::Vector3f(position.x(), position.y(),
			position.z());
	gadget::Vector3f b = field.getField(r / kpc);
	Vector3 bField = Vector3(b.x, b.y, b.z) * gauss;
	return bField;
}

void SPHMagneticField::updateSimulationVolume(const Vector3 &origin,
		double size) {
	gadget::Vector3f v = gadget::Vector3f(origin.x(), origin.y(), origin.z())
			/ kpc;
	field.init(v, size / kpc, database);
}

SPHMagneticFieldGrid::SPHMagneticFieldGrid(Vector3 origin, double size,
		size_t samples, const std::string filename) :
		field(samples) {
	gadget::Vector3f v = gadget::Vector3f(origin.x(), origin.y(), origin.z())
			/ kpc;
	field.init(v, size / kpc, database);
}

Vector3 SPHMagneticFieldGrid::getField(const Vector3 &position) const {
	gadget::Vector3f r = gadget::Vector3f(position.x(), position.y(),
			position.z());
	gadget::Vector3f b = field.getField(r / kpc);
	Vector3 bField = Vector3(b.x, b.y, b.z) * gauss;
	return bField;
}

void SPHMagneticFieldGrid::updateSimulationVolume(const Vector3 &origin,
		double size) {
	gadget::Vector3f v = gadget::Vector3f(origin.x(), origin.y(), origin.z())
			/ kpc;
	field.init(v, size / kpc, database);
}

} // namespace mpc


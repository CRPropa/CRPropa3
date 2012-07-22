#include "mpc/magneticField/SPHMagneticField.h"
#include "mpc/Units.h"

#include "kiss/logger.h"

namespace mpc {

SPHMagneticField::SPHMagneticField(Vector3d origin, double size, size_t samples,
		std::string filename) :
		field(samples) {
	database.open(filename);
	gadget::Vector3f v = gadget::Vector3f(origin.x, origin.y, origin.z) / kpc;
	field.init(v, size / kpc, database);
}

SPHMagneticField::SPHMagneticField(size_t samples, std::string filename) :
		field(samples) {
	database.open(filename);
}

Vector3d SPHMagneticField::getField(const Vector3d &position) const {
	gadget::Vector3f r = gadget::Vector3f(position.x, position.y, position.z);
	gadget::Vector3f b;
	bool isGood = field.getField(r / kpc, b);
	if (!isGood)
			std::cout << "mpc::SPHMagneticField invalid position : " << position << std::endl;
	Vector3d bField = Vector3d(b.x, b.y, b.z) * gauss;
	return bField;
}

double SPHMagneticField::getRho(const Vector3d& position) const {
	gadget::Vector3f r = gadget::Vector3f(position.x, position.y, position.z);
	size_t overlaps = 0;
	float rho;
	bool isGood = field.getRho(r / kpc, overlaps, rho);
	if (!isGood)
		std::cout << "mpc::SPHMagneticField invalid position : " << position << std::endl;
	return rho * 1.98892e40 * kilogram * pow(0.7, 2) / pow(kpc, 3);
}

SPHMagneticFieldGrid::SPHMagneticFieldGrid(Vector3d origin, double size,
		size_t samples, std::string filename) :
		samples(samples), field(samples), cacheEnabled(false) {
	database.open(filename);
	gadget::Vector3f v = gadget::Vector3f(origin.x, origin.y, origin.z) / kpc;
	field.init(v, size / kpc, database);
}

SPHMagneticFieldGrid::SPHMagneticFieldGrid(const size_t samples,
		const std::string filename) :
		samples(samples), field(samples) {
	database.open(filename);
}

Vector3d SPHMagneticFieldGrid::getField(const Vector3d &position) const {
	gadget::Vector3f r = gadget::Vector3f(position.x, position.y, position.z);
	gadget::Vector3f b;
	bool isGood = field.getField(r / kpc, b);
	if (!isGood)
			std::cout << "mpc::SPHMagneticField invalid position : " << position << std::endl;
	Vector3d bField = Vector3d(b.x, b.y, b.z) * gauss;
	return bField;
}

void SPHMagneticFieldGrid::setCachePrefix(std::string prefix) {
	cachePrefix = prefix;
}

void SPHMagneticFieldGrid::setCacheEnabled(bool enabled) {
	cacheEnabled = enabled;
}

} // namespace mpc


#include "mpc/magneticField/sphMagneticField.h"
#include "mpc/Units.h"

namespace mpc {

SPHMagneticField::SPHMagneticField(Vector3 origin, double size, size_t n) {
	gadget::Vector3f v = gadget::Vector3f(origin.x(), origin.y(), origin.z()) / kpc;
	gadgetField.reset(new gadget::DirectMagneticField(v, size / kpc));
	gadgetField->init(n);
}

Vector3 SPHMagneticField::getField(const Vector3 &position) const {
	gadget::Vector3f r = gadget::Vector3f(position.x(), position.y(), position.z());
	gadget::Vector3f b = gadgetField->getField(r / kpc);
	Vector3 bField = Vector3(b.x, b.y, b.z) * gauss;
	return bField;
}

SPHMagneticFieldGrid::SPHMagneticFieldGrid(Vector3 origin, double size, size_t n) {
	gadget::Vector3f v = gadget::Vector3f(origin.x(), origin.y(), origin.z()) / kpc;
	gadgetField.reset(new gadget::SampledMagneticField(v, size / kpc));
	gadgetField->init(size / n);
}

Vector3 SPHMagneticFieldGrid::getField(const Vector3 &position) const {
	gadget::Vector3f r = gadget::Vector3f(position.x(), position.y(), position.z());
	gadget::Vector3f b = gadgetField->getField(r / kpc);
	Vector3 bField = Vector3(b.x, b.y, b.z) * gauss;
	return bField;
}

} // namespace mpc


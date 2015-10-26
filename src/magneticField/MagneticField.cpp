#include "crpropa/magneticField/MagneticField.h"

namespace crpropa {

PeriodicMagneticField::PeriodicMagneticField(ref_ptr<MagneticField> field,
		const Vector3d &extends) :
		field(field), extends(extends), origin(0, 0, 0), reflective(false) {

}

PeriodicMagneticField::PeriodicMagneticField(ref_ptr<MagneticField> field,
		const Vector3d &extends, const Vector3d &origin, bool reflective) :
		field(field), extends(extends), origin(origin), reflective(reflective) {

}

Vector3d &PeriodicMagneticField::getOrigin() {
	return origin;
}

void PeriodicMagneticField::setOrigin(const Vector3d &origin) {
	this->origin = origin;
}

Vector3d &PeriodicMagneticField::getExtends() {
	return extends;
}

void PeriodicMagneticField::setExtends(const Vector3d &origin) {
	this->extends = extends;
}

bool PeriodicMagneticField::isReflective() {
	return reflective;
}

void PeriodicMagneticField::setReflective(bool reflective) {
	this->reflective = reflective;
}

Vector3d PeriodicMagneticField::getField(const Vector3d &position) const {
	Vector3d n = ((position - origin) / extends).floor();
	Vector3d p = position - origin - n * extends;

	if (reflective) {
		long mx = (long) ::fabs(n.x) % 2;
		if (mx == 1)
			p.x = extends.x - p.x;
		long my = (long) ::fabs(n.y) % 2;
		if (my == 1)
			p.y = extends.y - p.y;
		long mz = (long) ::fabs(n.z) % 2;
		if (mz == 1)
			p.z = extends.z - p.z;
	}

	return field->getField(p);
}

void MagneticFieldList::addField(ref_ptr<MagneticField> field) {
	fields.push_back(field);
}

Vector3d MagneticFieldList::getField(const Vector3d &position) const {
	Vector3d b;
	for (int i = 0; i < fields.size(); i++)
		b += fields[i]->getField(position);
	return b;
}

MagneticFieldEvolution::MagneticFieldEvolution(ref_ptr<MagneticField> field,
	double m) :
	field(field), m(m) {
}

Vector3d MagneticFieldEvolution::getField(const Vector3d &position,
	double z) const {
	return field->getField(position) * pow(1+z, m);
}

} // namespace crpropa

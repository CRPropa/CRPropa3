#include "mpc/magneticField/MagneticField.h"

namespace mpc {

void MagneticFieldList::addField(ref_ptr<MagneticField> field) {
	fields.push_back(field);
}

Vector3d MagneticFieldList::getField(const Vector3d &position) const {
	Vector3d b;
	for (int i = 0; i < fields.size(); i++)
		b += fields[i]->getField(position);
	return b;
}

} // namespace mpc

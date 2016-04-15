#ifndef CRPROPA_QUIMBYMAGNETICFIELD_H
#define CRPROPA_QUIMBYMAGNETICFIELD_H

#ifdef CRPROPA_HAVE_QUIMBY

#include "crpropa/Units.h"
#include "crpropa/magneticField/MagneticField.h"

#include "quimby/MagneticField.h"

namespace crpropa {

/**
 @class QuimbyMagneticField
 @brief Wrapper for quimby::MagneticField
 */
class QuimbyMagneticField: public MagneticField {
	quimby::ref_ptr<quimby::MagneticField> field;
public:
	QuimbyMagneticField(quimby::ref_ptr<quimby::MagneticField> field) : field(field) {

	}
	QuimbyMagneticField(quimby::MagneticField *field) : field(field) {
	}
	Vector3d getField(const Vector3d &position) const {
		quimby::Vector3f b, r = quimby::Vector3f(position.x, position.y, position.z);
		bool isGood = field->getField(r / kpc, b);
		if (!isGood)
			std::cerr << "crpropa::SPHMagneticField invalid position : " << position
					<< std::endl;
		return Vector3d(b.x, b.y, b.z) * gauss;
	}
};
#if 1
/**
 @class QuimbyMagneticFieldAdapter
 @brief Wrapper to use crpropa::MagneticField in Quimby
 */
class QuimbyMagneticFieldAdapter: public quimby::MagneticField {
	crpropa::ref_ptr<crpropa::MagneticField> field;
public:
	QuimbyMagneticFieldAdapter(crpropa::ref_ptr<crpropa::MagneticField> field) : field(field) {

	}

	bool getField(const quimby::Vector3f &position, quimby::Vector3f &b) const {
		crpropa::Vector3d r = crpropa::Vector3d(position.x, position.y, position.z) * crpropa::kpc;
		crpropa::Vector3d B = field->getField(r);
		b = quimby::Vector3f(B.x, B.y, B.z) / gauss;
		return true;
	}
};
#endif

} // namespace crpropa



#endif // CRPROPA_HAVE_QUIMBY
#endif // CRPROPA_QUIMBYMAGNETICFIELD_H

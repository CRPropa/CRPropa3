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
	Vector3d getField(const Vector3d &position) const {
		quimby::Vector3f b, r = quimby::Vector3f(position.x, position.y, position.z);
		bool isGood = field->getField(r / kpc, b);
		if (!isGood)
			std::cerr << "mpc::SPHMagneticField invalid position : " << position
					<< std::endl;
		return Vector3d(b.x, b.y, b.z) * gauss;
	}
};

} // namespace crpropaCRPROPA

#endif // CRPROPA_HAVE_QUIMBY
#endif // CRPROPA_QUIMBYMAGNETICFIELD_H

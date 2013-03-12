#ifndef MPC_GADGETMAGNETICFIELD_H_
#define MPC_GADGETMAGNETICFIELD_H_

#ifdef MPC_HAVE_GADGET

#include "mpc/Units.h"
#include "mpc/magneticField/MagneticField.h"

#include "gadget/MagneticField.h"

namespace mpc {

/**
 @class GadgetMagneticField
 @brief Wrapper for gadget::MagneticField
 */
class GadgetMagneticField: public MagneticField {
	gadget::ref_ptr<gadget::MagneticField> field;
public:
	GadgetMagneticField(gadget::ref_ptr<gadget::MagneticField> field) : field(field) {

	}
	Vector3d getField(const Vector3d &position) const {
		gadget::Vector3f b, r = gadget::Vector3f(position.x, position.y, position.z);
		bool isGood = field->getField(r / kpc, b);
		if (!isGood)
			std::cerr << "mpc::SPHMagneticField invalid position : " << position
					<< std::endl;
		return Vector3d(b.x, b.y, b.z) * gauss;
	}
};

} // namespace mpc

#endif // MPC_HAVE_GADGET
#endif // MPC_GADGETMAGNETICFIELD_H_

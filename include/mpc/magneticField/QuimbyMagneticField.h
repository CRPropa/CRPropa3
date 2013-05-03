#ifndef MPC_QUIMBYMAGNETICFIELD_H_
#define MPC_QUIMBYMAGNETICFIELD_H_

#ifdef MPC_HAVE_QUIMBY

#include "mpc/Units.h"
#include "mpc/magneticField/MagneticField.h"

#include "quimby/MagneticField.h"

namespace mpc {

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

} // namespace mpc

#endif // MPC_HAVE_QUIMBY
#endif // MPC_QUIMBYMAGNETICFIELD_H_

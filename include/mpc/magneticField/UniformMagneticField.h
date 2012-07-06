#ifndef MPC_UNIFORMAGNETICFIELD_H_
#define MPC_UNIFORMAGNETICFIELD_H_

#include "mpc/magneticField/MagneticField.h"

namespace mpc {

/**
 @class UniformMagneticField
 @brief Magnetic field with one B-field vector.
 */
class UniformMagneticField: public MagneticField {
public:
	UniformMagneticField(const Vector3d &value) :
			value(value) {
	}
	Vector3d getField(const Vector3d &position) const {
		return value;
	}

	void updateSimulationVolume(const Vector3d &origin, double size) {

	}

private:
	Vector3d value;
};

} // namespace mpc

#endif /* MPC_UNIFORMAGNETICFIELD_H_ */

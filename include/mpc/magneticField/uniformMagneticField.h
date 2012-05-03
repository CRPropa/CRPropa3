#ifndef UNIFORMAGNETICFIELD_H_
#define UNIFORMAGNETICFIELD_H_

#include "mpc/magneticField/magneticField.h"
#include "mpc/Vector3.h"

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

#endif /* UNIFORMAGNETICFIELD_H_ */

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
	UniformMagneticField(const Vector3 &value) :
			value(value) {
	}
	Vector3 getField(const Vector3 &position) const {
		return value;
	}

private:
	Vector3 value;
};

} // namespace mpc

#endif /* UNIFORMAGNETICFIELD_H_ */

#ifndef UNIFORMAGNETICFIELD_HPP_
#define UNIFORMAGNETICFIELD_HPP_

#include "mpc/magneticField/magneticField.hpp"
#include "mpc/Vector3.h"

namespace mpc {

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

#endif /* UNIFORMAGNETICFIELD_HPP_ */

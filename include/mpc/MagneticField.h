#ifndef MAGNETICFIELD_H_
#define MAGNETICFIELD_H_

#include "mpc/Vector3.h"
#include "mpc/Units.h"

namespace mpc {

class MagneticField {
public:
	virtual ~MagneticField() {
	}
	virtual Vector3 getField(const Vector3 &position) const = 0;
};

class HomogeneousMagneticField: public MagneticField {
public:
	HomogeneousMagneticField(const Vector3 &value) :
			value(value) {
	}
	Vector3 getField(const Vector3 &position) const {
		return value;
	}
private:
	Vector3 value;
};

} // namespace mpc

#endif /* MAGNETICFIELD_H_ */

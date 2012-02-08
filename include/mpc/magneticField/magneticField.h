#ifndef MAGNETICFIELD_H_
#define MAGNETICFIELD_H_

#include "mpc/Vector3.h"

namespace mpc {

/**
	@class MagneticField
	@brief Magnetic field base class.
*/
class MagneticField {
public:
	virtual ~MagneticField() {
	}
	virtual Vector3 getField(const Vector3 &position) const = 0;
};

} // namespace mpc

#endif /* MAGNETICFIELD_H_ */

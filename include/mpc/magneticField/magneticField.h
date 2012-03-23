#ifndef MAGNETICFIELD_H_
#define MAGNETICFIELD_H_

#include "mpc/Vector3.h"
#include "mpc/Referenced.h"

namespace mpc {

/**
 @class MagneticField
 @brief Abstract base class for magnetic fields.
 */
class MagneticField: public Referenced {
public:
	virtual ~MagneticField() {
	}
	virtual Vector3 getField(const Vector3 &position) const = 0;
	virtual void updateSimulationVolume(const Vector3 &origin, double size) = 0;
};

} // namespace mpc

#endif /* MAGNETICFIELD_H_ */

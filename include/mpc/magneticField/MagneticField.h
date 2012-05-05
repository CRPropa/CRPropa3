#ifndef MPC_MAGNETICFIELD_H_
#define MPC_MAGNETICFIELD_H_

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
	virtual Vector3d getField(const Vector3d &position) const = 0;
	virtual void updateSimulationVolume(const Vector3d &origin, double size) = 0;
};

} // namespace mpc

#endif /* MPC_MAGNETICFIELD_H_ */

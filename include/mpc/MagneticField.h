#ifndef MAGNETICFIELD_H_
#define MAGNETICFIELD_H_

#include "mpc/ThreeVector.h"
#include "mpc/Units.h"

namespace mpc{

class MagneticField {
public:
	virtual Hep3Vector getField(const Hep3Vector &position) const = 0;
};

class HomogeneousMagneticField: public MagneticField {
public:
	HomogeneousMagneticField(const Hep3Vector &value) :
			value(value) {
	}
	Hep3Vector getField(const Hep3Vector &position) const {
		return value;
	}
private:
	Hep3Vector value;
};

} // namespace mpc

#endif /* MAGNETICFIELD_H_ */

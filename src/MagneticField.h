
#ifndef MAGNETICFIELD_H_
#define MAGNETICFIELD_H_

#include "ThreeVector.h"

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

#endif /* MAGNETICFIELD_H_ */

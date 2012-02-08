#ifndef _MAGNETICFIELDRING_H_
#define _MAGNETICFIELDRING_H_

#include "mpc/magneticField/magneticField.h"
#include <math.h>

namespace mpc {

class MagneticFieldRing: public MagneticField {

protected:
	Vector3 center;
	double Bnorm;
	double innerEdge;
	double outerEdge;
	double scaleHeight;

public:
	MagneticFieldRing(Vector3 center, double Bnorm, double innerEdge,
			double outerEdge, double scaleHeight) {
		this->center = center;
		this->Bnorm = Bnorm;
		this->innerEdge = innerEdge;
		this->outerEdge = outerEdge;
		this->scaleHeight = scaleHeight;
	}

	~MagneticFieldRing() {
	}

	Vector3 getField(const Vector3 &position) const {
		Vector3 b(0, 0, 0);
		Vector3 r = position - center;

		if ((r.mag() < innerEdge) || (r.mag() > outerEdge))
			return b;

		double phi = -1. * atan2(r.y(), r.x());
		double weight = exp(-fabs(r.z()) / scaleHeight);
		b.setX(Bnorm * sin(phi) * weight);
		b.setY(Bnorm * cos(phi) * weight);
		return b;
	}
};

} // namespace mpc

#endif // _MAGNETICFIELDRING_H_

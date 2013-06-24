#ifndef PHASEPOINT_H
#define PHASEPOINT_H

#include "crpropa/Vector3.h"

namespace crpropa {

/**
 @class PhasePoint
 @brief 6D point in (position - momentum) phase space
 */
class PhasePoint {
public:
	Vector3d a, b;

	PhasePoint() {
	}

	PhasePoint(const Vector3d &a, const Vector3d &b) :
			a(a), b(b) {
	}

	PhasePoint(double f) :
			a(Vector3d(f, f, f)), b(Vector3d(f, f, f)) {
	}

	PhasePoint operator *(double f) const {
		return PhasePoint(a * f, b * f);
	}

	PhasePoint &operator =(double f) {
		a = Vector3d(f, f, f);
		b = Vector3d(f, f, f);
		return *this;
	}

	PhasePoint &operator +=(const PhasePoint &p) {
		a += p.a;
		b += p.b;
		return *this;
	}

	PhasePoint operator +(const PhasePoint &p) const {
		return PhasePoint(a + p.a, b + p.b);
	}

	PhasePoint abs() const {
		return PhasePoint(a.abs(), b.abs());
	}
};

} // namespace crpropa

#endif /* PHASEPOINT_H */

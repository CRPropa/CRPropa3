#ifndef PHASEPOINT_H_
#define PHASEPOINT_H_

#include "mpc/Vector3.h"

namespace mpc {

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
		Vector3d abs_a( fabs(a.x), fabs(a.y), fabs(a.z) );
		Vector3d abs_b( fabs(b.x), fabs(b.y), fabs(b.z) );
		return PhasePoint(abs_a, abs_b);
	}
};

} // namespace mpc

#endif /* PHASEPOINT_H_ */

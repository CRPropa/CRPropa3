#ifndef PHASEPOINT_H_
#define PHASEPOINT_H_

#include "mpc/Vector3.h"

namespace mpc {

class PhasePoint {
public:
	Vector3 a, b;

	PhasePoint() {
	}

	PhasePoint(const Vector3 &a, const Vector3 &b) :
			a(a), b(b) {
	}

	PhasePoint(double f) :
			a(Vector3(f, f, f)), b(Vector3(f, f, f)) {
	}

	PhasePoint operator *(double f) const {
		return PhasePoint(a * f, b * f);
	}

	PhasePoint &operator =(double f) {
		a.set(f, f, f);
		b.set(f, f, f);
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
		Vector3 abs_a( fabs(a.x()), fabs(a.y()), fabs(a.z()) );
		Vector3 abs_b( fabs(b.x()), fabs(b.y()), fabs(b.z()) );
		return PhasePoint(abs_a, abs_b);
	}
};

} // namespace mpc

#endif /* PHASEPOINT_H_ */

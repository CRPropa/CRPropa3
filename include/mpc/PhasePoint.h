#ifndef PHASEPOINT_H_
#define PHASEPOINT_H_

#include "mpc/ThreeVector.h"

#include <limits>

namespace mpc {

class PhasePoint {
public:
	Hep3Vector a, b;

	PhasePoint() {
	}

	PhasePoint(const Hep3Vector &a, const Hep3Vector &b) :
			a(a), b(b) {
	}

	PhasePoint(double f) :
			a(Hep3Vector(f, f, f)), b(Hep3Vector(f, f, f)) {
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
		Hep3Vector abs_a( fabs(a.x()), fabs(a.y()), fabs(a.z()) );
		Hep3Vector abs_b( fabs(b.x()), fabs(b.y()), fabs(b.z()) );
		return PhasePoint(abs_a, abs_b);
	}
};
} // namespace mpc

#endif /* PHASEPOINT_H_ */

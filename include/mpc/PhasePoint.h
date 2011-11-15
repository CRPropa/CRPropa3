#ifndef PHASEPOINT_H_
#define PHASEPOINT_H_

#include "ThreeVector.h"
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
		Hep3Vector abs_a(fabs(a.x()), fabs(a.y()), fabs(a.z()));
		Hep3Vector abs_b(fabs(b.x()), fabs(b.y()), fabs(b.z()));
		return PhasePoint(abs_a, abs_b);
	}

	double getMaxRatio(const PhasePoint &rel) {
		double m = -std::numeric_limits<double>::infinity();
		m = std::max(m, fabs(a.x() / rel.a.x()));
		m = std::max(m, fabs(a.y() / rel.a.y()));
		m = std::max(m, fabs(a.z() / rel.a.z()));
		m = std::max(m, fabs(b.x() / rel.b.x()));
		m = std::max(m, fabs(b.y() / rel.b.y()));
		m = std::max(m, fabs(b.z() / rel.b.z()));
		return m;
	}
};
} // namespace mpc

#endif /* PHASEPOINT_H_ */

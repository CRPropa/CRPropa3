// -*- C++ -*-
// $Id: ThreeVector.cc,v 1.3 2003/08/13 20:00:14 garren Exp $
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of the Vector3 class.
//
// See also ThreeVectorR.cc for implementation of Vector3 methods which
// would couple in all the HepRotation methods.
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "mpc/Vector3.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace mpc {
void Vector3::setMag(double ma) {
	double factor = mag();
	if (factor == 0) {
		std::cerr << "Vector3::setMag : zero vector can't be stretched"
				<< std::endl;
	} else {
		factor = ma / factor;
		setX(x() * factor);
		setY(y() * factor);
		setZ(z() * factor);
	}
}

double Vector3::operator ()(int i) const {
	switch (i) {
	case X:
		return x();
	case Y:
		return y();
	case Z:
		return z();
	default:
		std::cerr << "Vector3 subscripting: bad index (" << i << ")"
				<< std::endl;
	}
	return 0.;
}

double & Vector3::operator ()(int i) {
	static double dummy;
	switch (i) {
	case X:
		return dx;
	case Y:
		return dy;
	case Z:
		return dz;
	default:
		std::cerr << "Vector3 subscripting: bad index (" << i << ")"
				<< std::endl;
		return dummy;
	}
}

Vector3 & Vector3::rotateUz(const Vector3& NewUzVector) {
	// NewUzVector must be normalized !

	double u1 = NewUzVector.x();
	double u2 = NewUzVector.y();
	double u3 = NewUzVector.z();
	double up = u1 * u1 + u2 * u2;

	if (up > 0) {
		up = sqrt(up);
		double px = dx, py = dy, pz = dz;
		dx = (u1 * u3 * px - u2 * py) / up + u1 * pz;
		dy = (u2 * u3 * px + u1 * py) / up + u2 * pz;
		dz = -up * px + u3 * pz;
	} else if (u3 < 0.) {
		dx = -dx;
		dz = -dz;
	} // phi=0  teta=pi
	else {
	};
	return *this;
}

double Vector3::pseudoRapidity() const {
	double m = mag();
	if (m == 0)
		return 0.0;
	if (m == z())
		return 1.0E72;
	if (m == -z())
		return -1.0E72;
	return 0.5 * log((m + z()) / (m - z()));
}

std::ostream & operator<<(std::ostream & os, const Vector3 & v) {
	return os << "(" << v.x() << "," << v.y() << "," << v.z() << ")";
}

std::istream & operator>>(std::istream & is, Vector3 & v) {
	double x, y, z;
	is >> x;
	is >> y;
	is >> z;
	v.set(x, y, z);
	return is;
} // operator>>()

const Vector3 HepXHat(1.0, 0.0, 0.0);
const Vector3 HepYHat(0.0, 1.0, 0.0);
const Vector3 HepZHat(0.0, 0.0, 1.0);

//-------------------
//
// New methods introduced when ZOOM PhysicsVectors was merged in:
//
//-------------------

Vector3 & Vector3::rotateX(double phi) {
	double sinphi = sin(phi);
	double cosphi = cos(phi);
	double ty;
	ty = dy * cosphi - dz * sinphi;
	dz = dz * cosphi + dy * sinphi;
	dy = ty;
	return *this;
} /* rotateX */

Vector3 & Vector3::rotateY(double phi) {
	double sinphi = sin(phi);
	double cosphi = cos(phi);
	double tz;
	tz = dz * cosphi - dx * sinphi;
	dx = dx * cosphi + dz * sinphi;
	dz = tz;
	return *this;
} /* rotateY */

Vector3 & Vector3::rotateZ(double phi) {
	double sinphi = sin(phi);
	double cosphi = cos(phi);
	double tx;
	tx = dx * cosphi - dy * sinphi;
	dy = dy * cosphi + dx * sinphi;
	dx = tx;
	return *this;
} /* rotateZ */

bool Vector3::isNear(const Vector3 & v, double epsilon) const {
	double limit = dot(v) * epsilon * epsilon;
	return ((*this - v).mag2() <= limit);
} /* isNear() */

double Vector3::howNear(const Vector3 & v) const {
	// | V1 - V2 | **2  / V1 dot V2, up to 1
	double d = (*this - v).mag2();
	double vdv = dot(v);
	if ((vdv > 0) && (d < vdv)) {
		return sqrt(d / vdv);
	} else if ((vdv == 0) && (d == 0)) {
		return 0;
	} else {
		return 1;
	}
} /* howNear */

double Vector3::deltaPhi(const Vector3 & v2) const {
	double dphi = v2.getPhi() - getPhi();
	if (dphi > M_PI) {
		dphi -= 2 * M_PI;
	} else if (dphi <= -M_PI) {
		dphi += 2 * M_PI;
	}
	return dphi;
} /* deltaPhi */

double Vector3::deltaR(const Vector3 & v) const {
	double a = eta() - v.eta();
	double b = deltaPhi(v);
	return sqrt(a * a + b * b);
} /* deltaR */

double Vector3::cosTheta(const Vector3 & q) const {
	double arg;
	double ptot2 = mag2() * q.mag2();
	if (ptot2 <= 0) {
		arg = 0.0;
	} else {
		arg = dot(q) / sqrt(ptot2);
		if (arg > 1.0)
			arg = 1.0;
		if (arg < -1.0)
			arg = -1.0;
	}
	return arg;
}

double Vector3::cos2Theta(const Vector3 & q) const {
	double arg;
	double ptot2 = mag2();
	double qtot2 = q.mag2();
	if (ptot2 == 0 || qtot2 == 0) {
		arg = 1.0;
	} else {
		double pdq = dot(q);
		arg = (pdq / ptot2) * (pdq / qtot2);
		// More naive methods overflow on vectors which can be squared
		// but can't be raised to the 4th power.
		if (arg > 1.0)
			arg = 1.0;
	}
	return arg;
}

void Vector3::setEta(double eta) {
	double phi = 0;
	double r;
	if ((dx == 0) && (dy == 0)) {
		if (dz == 0) {
			std::cerr
					<< "Attempt to set eta of zero vector -- vector is unchanged"
					<< std::endl;
			return;
		}
		std::cerr
				<< "Attempt to set eta of vector along Z axis -- will use phi = 0"
				<< std::endl;
		r = fabs(dz);
	} else {
		r = getR();
		phi = getPhi();
	}
	double tanHalfTheta = exp(-eta);
	double cosTheta = (1 - tanHalfTheta * tanHalfTheta)
			/ (1 + tanHalfTheta * tanHalfTheta);
	dz = r * cosTheta;
	double rho = r * sqrt(1 - cosTheta * cosTheta);
	dy = rho * sin(phi);
	dx = rho * cos(phi);
	return;
}

void Vector3::setCylTheta(double theta) {

	// In cylindrical coords, set theta while keeping rho and phi fixed

	if ((dx == 0) && (dy == 0)) {
		if (dz == 0) {
			std::cerr
					<< "Attempt to set cylTheta of zero vector -- vector is unchanged"
					<< std::endl;
			return;
		}
		if (theta == 0) {
			dz = fabs(dz);
			return;
		}
		if (theta == M_PI) {
			dz = -fabs(dz);
			return;
		}
		std::cerr << "Attempt set cylindrical theta of vector along Z axis "
				"to a non-trivial value, while keeping rho fixed -- "
				"will return zero vector" << std::endl;
		dz = 0;
		return;
	}
	if ((theta < 0) || (theta > M_PI)) {
		std::cerr
				<< "Setting Cyl theta of a vector based on a value not in [0, PI]"
				<< std::endl;
		// No special return needed if warning is ignored.
	}
	double phi(getPhi());
	double rho = getRho();
	if ((theta == 0) || (theta == M_PI)) {
		std::cerr << "Attempt to set cylindrical theta to 0 or PI "
				"while keeping rho fixed -- infinite Z will be computed"
				<< std::endl;
		dz = (theta == 0) ? 1.0E72 : -1.0E72;
		return;
	}
	dz = rho / tan(theta);
	dy = rho * sin(phi);
	dx = rho * cos(phi);

} /* setCylTheta */

void Vector3::setCylEta(double eta) {

	// In cylindrical coords, set eta while keeping rho and phi fixed

	double theta = 2 * atan(exp(-eta));

	//-| The remaining code is similar to setCylTheta,  The reason for
	//-| using a copy is so as to be able to change the messages in the
	//-| ZMthrows to say eta rather than theta.  Besides, we assumedly
	//-| need not check for theta of 0 or PI.

	if ((dx == 0) && (dy == 0)) {
		if (dz == 0) {
			std::cerr
					<< "Attempt to set cylEta of zero vector -- vector is unchanged"
					<< std::endl;
			return;
		}
		if (theta == 0) {
			dz = fabs(dz);
			return;
		}
		if (theta == M_PI) {
			dz = -fabs(dz);
			return;
		}
		std::cerr << "Attempt set cylindrical eta of vector along Z axis "
				"to a non-trivial value, while keeping rho fixed -- "
				"will return zero vector" << std::endl;
		dz = 0;
		return;
	}
	double phi(getPhi());
	double rho = getRho();
	dz = rho / tan(theta);
	dy = rho * sin(phi);
	dx = rho * cos(phi);

} /* setCylEta */

Vector3 operator/(const Vector3 & v1, double c) {
	if (c == 0) {
		std::cerr << "Attempt to divide vector by 0 -- "
				"will produce infinities and/or NANs" << std::endl;
	}
	double oneOverC = 1.0 / c;
	return Vector3(v1.x() * oneOverC, v1.y() * oneOverC, v1.z() * oneOverC);
} /* v / c */

Vector3 & Vector3::operator/=(double c) {
	if (c == 0) {
		std::cerr << "Attempt to do vector /= 0 -- "
				"division by zero would produce infinite or NAN components"
				<< std::endl;
	}
	double oneOverC = 1.0 / c;
	dx *= oneOverC;
	dy *= oneOverC;
	dz *= oneOverC;
	return *this;
}

void Vector3::setSpherical(double r, double theta, double phi) {
	if (r < 0) {
		throw std::runtime_error("Spherical coordinates set with negative   R");
		// No special return needed if warning is ignored.
	}
	if ((theta < 0) || (theta > M_PI)) {
		throw std::runtime_error(
				"Spherical coordinates set with theta not in [0, PI]");
		// No special return needed if warning is ignored.
	}
	dz = r * cos(theta);
	double rho(r * sin(theta));
	dy = rho * sin(phi);
	dx = rho * cos(phi);
	return;
} /* setSpherical (r, theta, phi) */

void Vector3::setCylindrical(double rho, double phi, double z) {
	if (rho < 0) {
		throw std::runtime_error(
				"Cylindrical coordinates supplied with negative Rho");
		// No special return needed if warning is ignored.
	}
	dz = z;
	dy = rho * sin(phi);
	dx = rho * cos(phi);
	return;
} /* setCylindrical (r, phi, z) */

void Vector3::setRhoPhiTheta(double rho, double phi, double theta) {
	if (rho == 0) {
		throw std::runtime_error(
				"Attempt set vector components rho, phi, theta with zero rho -- "
						"zero vector is returned, ignoring theta and phi");
		dx = 0;
		dy = 0;
		dz = 0;
		return;
	}
	if ((theta == 0) || (theta == M_PI)) {
		throw std::runtime_error(
				"Attempt set cylindrical vector vector with finite rho and "
						"theta along the Z axis:  infinite Z would be computed");
	}
	if ((theta < 0) || (theta > M_PI)) {
		throw std::runtime_error(
				"Rho, phi, theta set with theta not in [0, PI]");
		// No special return needed if warning is ignored.
	}
	dz = rho / tan(theta);
	dy = rho * sin(phi);
	dx = rho * cos(phi);
	return;
} /* setCyl (rho, phi, theta) */

void Vector3::setRhoPhiEta(double rho, double phi, double eta) {
	if (rho == 0) {
		throw std::runtime_error(
				"Attempt set vector components rho, phi, eta with zero rho -- "
						"zero vector is returned, ignoring eta and phi");
		dx = 0;
		dy = 0;
		dz = 0;
		return;
	}
	double theta(2 * atan(exp(-eta)));
	dz = rho / tan(theta);
	dy = rho * sin(phi);
	dx = rho * cos(phi);
	return;
} /* setCyl (rho, phi, eta) */

//************
//    - 3 -
// Comparisons
//
//************

int Vector3::compare(const Vector3 & v) const {
	if (dz > v.dz) {
		return 1;
	} else if (dz < v.dz) {
		return -1;
	} else if (dy > v.dy) {
		return 1;
	} else if (dy < v.dy) {
		return -1;
	} else if (dx > v.dx) {
		return 1;
	} else if (dx < v.dx) {
		return -1;
	} else {
		return 0;
	}
} /* Compare */

bool Vector3::operator >(const Vector3 & v) const {
	return (compare(v) > 0);
}
bool Vector3::operator <(const Vector3 & v) const {
	return (compare(v) < 0);
}
bool Vector3::operator>=(const Vector3 & v) const {
	return (compare(v) >= 0);
}
bool Vector3::operator<=(const Vector3 & v) const {
	return (compare(v) <= 0);
}

//-********
// Nearness
//-********

// These methods all assume you can safely take mag2() of each vector.
// Absolutely safe but slower and much uglier alternatives were
// provided as build-time options in ZOOM SpaceVectors.
// Also, much smaller codes were provided tht assume you can square
// mag2() of each vector; but those return bad answers without warning
// when components exceed 10**90.
//
// IsNear, HowNear, and DeltaR are found in ThreeVector.cc

double Vector3::howParallel(const Vector3 & v) const {
	// | V1 x V2 | / | V1 dot V2 |
	double v1v2 = fabs(dot(v));
	if (v1v2 == 0) {
		// Zero is parallel to no other vector except for zero.
		return ((mag2() == 0) && (v.mag2() == 0)) ? 0 : 1;
	}
	Vector3 v1Xv2(cross(v));
	double abscross = v1Xv2.mag();
	if (abscross >= v1v2) {
		return 1;
	} else {
		return abscross / v1v2;
	}
} /* howParallel() */

bool Vector3::isParallel(const Vector3 & v, double epsilon) const {
	// | V1 x V2 | **2  <= epsilon **2 | V1 dot V2 | **2
	// V1 is *this, V2 is v

	static const double TOOBIG = pow(2.0, 507);
	static const double SCALE = pow(2.0, -507);
	double v1v2 = fabs(dot(v));
	if (v1v2 == 0) {
		return ((mag2() == 0) && (v.mag2() == 0));
	}
	if (v1v2 >= TOOBIG) {
		Vector3 sv1(*this * SCALE);
		Vector3 sv2(v * SCALE);
		Vector3 sv1Xsv2 = sv1.cross(sv2);
		double x2 = sv1Xsv2.mag2();
		double limit = v1v2 * SCALE * SCALE;
		limit = epsilon * epsilon * limit * limit;
		return (x2 <= limit);
	}

	// At this point we know v1v2 can be squared.

	Vector3 v1Xv2(cross(v));
	if ((fabs(v1Xv2.dx) > TOOBIG) || (fabs(v1Xv2.dy) > TOOBIG)
			|| (fabs(v1Xv2.dz) > TOOBIG)) {
		return false;
	}

	return ((v1Xv2.mag2()) <= ((epsilon * v1v2) * (epsilon * v1v2)));

} /* isParallel() */

double Vector3::howOrthogonal(const Vector3 & v) const {
	// | V1 dot V2 | / | V1 x V2 |

	double v1v2 = fabs(dot(v));
	//-| Safe because both v1 and v2 can be squared
	if (v1v2 == 0) {
		return 0; // Even if one or both are 0, they are considered orthogonal
	}
	Vector3 v1Xv2(cross(v));
	double abscross = v1Xv2.mag();
	if (v1v2 >= abscross) {
		return 1;
	} else {
		return v1v2 / abscross;
	}

} /* howOrthogonal() */

bool Vector3::isOrthogonal(const Vector3 & v, double epsilon) const {
// | V1 x V2 | **2  <= epsilon **2 | V1 dot V2 | **2
// V1 is *this, V2 is v

	static const double TOOBIG = pow(2.0, 507);
	static const double SCALE = pow(2.0, -507);
	double v1v2 = fabs(dot(v));
	//-| Safe because both v1 and v2 can be squared
	if (v1v2 >= TOOBIG) {
		Vector3 sv1(*this * SCALE);
		Vector3 sv2(v * SCALE);
		Vector3 sv1Xsv2 = sv1.cross(sv2);
		double x2 = sv1Xsv2.mag2();
		double limit = epsilon * epsilon * x2;
		double y2 = v1v2 * SCALE * SCALE;
		return (y2 * y2 <= limit);
	}

	// At this point we know v1v2 can be squared.

	Vector3 eps_v1Xv2(cross(epsilon * v));
	if ((fabs(eps_v1Xv2.dx) > TOOBIG) || (fabs(eps_v1Xv2.dy) > TOOBIG)
			|| (fabs(eps_v1Xv2.dz) > TOOBIG)) {
		return true;
	}

	// At this point we know all the math we need can be done.

	return (v1v2 * v1v2 <= eps_v1Xv2.mag2());

} /* isOrthogonal() */

double Vector3::setTolerance(double tol) {
// Set the tolerance for Vector3s to be considered near one another
	double oldTolerance(tolerance);
	tolerance = tol;
	return oldTolerance;
}

//-***********************
// Helper Methods:
//	negativeInfinity()
//-***********************

double Vector3::negativeInfinity() const {
	// A byte-order-independent way to return -Infinity
	struct Dib {
		union {
			double d;
			unsigned char i[8];
		} u;
	};
	Dib negOne;
	Dib posTwo;
	negOne.u.d = -1.0;
	posTwo.u.d = 2.0;
	Dib value;
	int k;
	for (k = 0; k < 8; k++) {
		value.u.i[k] = negOne.u.i[k] | posTwo.u.i[k];
	}
	return value.u.d;
}
double Vector3::beta() const {
	double b = sqrt(mag2());
	if (b >= 1) {
		throw std::runtime_error(
				"Beta taken for Vector3 of at least unit length");
	}
	return b;
}

double Vector3::gamma() const {
	double beta = sqrt(mag2());
	if (beta == 1) {
		throw std::runtime_error(
				"Gamma taken for Vector3 of unit magnitude -- infinite result");
	}
	if (beta > 1) {
		throw std::runtime_error(
				"Gamma taken for Vector3 of more than unit magnitude -- "
						"the sqrt function would return NAN");
	}
	return 1 / sqrt(1 - beta * beta);
}

double Vector3::rapidity() const {
	if (fabs(dz) == 1) {
		throw std::runtime_error(
				"Rapidity in Z direction taken for Vector3 with |Z| = 1 -- \n"
						"the log should return infinity");
	}
	if (fabs(dz) > 1) {
		throw std::runtime_error(
				"Rapidity in Z direction taken for Vector3 with |Z| > 1 -- \n"
						"the log would return a NAN");
	}
// Want inverse tanh(dz):
	return (.5 * log((1 + dz) / (1 - dz)));
}

double Vector3::coLinearRapidity() const {
	double b = beta();
	if (b == 1) {
		throw std::runtime_error(
				"Co-linear Rapidity taken for Vector3 of unit length -- "
						"the log should return infinity");
	}
	if (b > 1) {
		throw std::runtime_error(
				"Co-linear Rapidity taken for Vector3 of more than unit length -- "
						"the log would return a NAN");
	}
// Want inverse tanh(b):
	return (.5 * log((1 + b) / (1 - b)));
}

//-***********************************************
// Other properties relative to a reference vector
//-***********************************************

Vector3 Vector3::project(const Vector3 & v2) const {
	double mag2v2 = v2.mag2();
	if (mag2v2 == 0) {
		throw std::runtime_error(
				"Attempt to take projection of vector against zero reference vector ");
		return project();
	}
	return (v2 * (dot(v2) / mag2v2));
}

double Vector3::rapidity(const Vector3 & v2) const {
	double vmag = v2.mag();
	if (vmag == 0) {
		throw std::runtime_error("Rapidity taken with respect to zero vector");
		return 0;
	}
	double z = dot(v2) / vmag;
	if (fabs(z) >= 1) {
		throw std::runtime_error("Rapidity taken for too large a Vector3 "
				"-- would return infinity or NAN");
	}
// Want inverse tanh(z):
	return (.5 * log((1 + z) / (1 - z)));
}

double Vector3::eta(const Vector3 & v2) const {
	// Defined as    -log ( tan ( .5* theta(u) ) );
	//
	// Quicker is to use cosTheta:
	// tan (theta/2) = sin(theta)/(1 + cos(theta))

	double r = getR();
	double v2r = v2.mag();
	if ((r == 0) || (v2r == 0)) {
		throw std::runtime_error(
				"Cannot find pseudorapidity of a zero vector relative to a vector");
		return 0.;
	}
	double c = dot(v2) / (r * v2r);
	if (c >= 1) {
		c = 1; //-| We don't want to return NAN because of roundoff
		throw std::runtime_error(
				"Pseudorapidity of vector relative to parallel vector -- "
						"will give infinite result");
		// We can just go on; tangent will be 0, so
		// log (tangent) will be -INFINITY, so result
		// will be +INFINITY.
	}
	if (c <= -1) {
		throw std::runtime_error(
				"Pseudorapidity of vector relative to anti-parallel vector -- "
						"will give negative infinite result");
		//-| We don't want to return NAN because of roundoff
		return (negativeInfinity());
		//  If we just went on, the tangent would be NAN
		//  so return would be NAN.  But the proper limit
		// of tan is +Infinity, so the return should be
		// -INFINITY.
	}

	double tangent = sqrt(1 - c * c) / (1 + c);
	return (-log(tangent));

} /* eta (u) */

//-*********************************************
//			- 6 -
// Decomposition of an angle between two vectors
//
//-*********************************************

double Vector3::polarAngle(const Vector3 & v2) const {
	return fabs(v2.getTheta() - getTheta());
} /* polarAngle */

double Vector3::polarAngle(const Vector3 & v2, const Vector3 & ref) const {
	return fabs(v2.angle(ref) - angle(ref));
} /* polarAngle (v2, ref) */

// double Vector3::azimAngle (const Vector3 & v2) const
// is now in the .icc file as deltaPhi(v2)

double Vector3::azimAngle(const Vector3 & v2, const Vector3 & ref) const {

	Vector3 vperp(perpPart(ref));
	if (vperp.mag2() == 0) {
		throw std::runtime_error(
				"Cannot find azimuthal angle with reference direction parallel to "
						"vector 1 -- will return zero");
		return 0;
	}

	Vector3 v2perp(v2.perpPart(ref));
	if (v2perp.mag2() == 0) {
		throw std::runtime_error(
				"Cannot find azimuthal angle with reference direction parallel to "
						"vector 2 -- will return zero");
		return 0;
	}

	double ang = vperp.angle(v2perp);

	// Now compute the sign of the answer:  that of U*(VxV2) or
	// the equivalent expression V*(V2xU).

	if (dot(v2.cross(ref)) >= 0) {
		return ang;
	} else {
		return -ang;
	}

	//-| Note that if V*(V2xU) is zero, we want to return 0 or PI
	//-| depending on whether vperp is aligned or antialigned with v2perp.
	//-| The computed angle() expression does this properly.

} /* azimAngle (v2, ref) */
//-************************
// rotate about axis
//-************************

Vector3 & Vector3::rotate(const Vector3 & axis, double delta) {
	double r = axis.mag();
	if (r == 0) {
		throw std::runtime_error(
				"Attempt to rotate around a zero vector axis! ");
		return *this;
	}
	register double scale = 1.0 / r;
	register double ux = scale * axis.getX();
	register double uy = scale * axis.getY();
	register double uz = scale * axis.getZ();
	double cd = cos(delta);
	double sd = sin(delta);
	register double ocd = 1 - cd;
	double rx;
	double ry;
	double rz;

	{
		register double ocdux = ocd * ux;
		rx = dx * (cd + ocdux * ux) + dy * (ocdux * uy - sd * uz)
				+ dz * (ocdux * uz + sd * uy);
	}

	{
		register double ocduy = ocd * uy;
		ry = dy * (cd + ocduy * uy) + dz * (ocduy * uz - sd * ux)
				+ dx * (ocduy * ux + sd * uz);
	}

	{
		register double ocduz = ocd * uz;
		rz = dz * (cd + ocduz * uz) + dx * (ocduz * ux - sd * uy)
				+ dy * (ocduz * uy + sd * ux);
	}

	dx = rx;
	dy = ry;
	dz = rz;

	return *this;
} /* rotate */

//-****************************
// rotate by three euler angles
//-****************************

Vector3 & Vector3::rotate(double phi, double theta, double psi) {

	double rx;
	double ry;
	double rz;

	register double sinPhi = sin(phi), cosPhi = cos(phi);
	register double sinTheta = sin(theta), cosTheta = cos(theta);
	register double sinPsi = sin(psi), cosPsi = cos(psi);

	rx = (cosPsi * cosPhi - cosTheta * sinPsi * sinPhi) * dx
			+ (cosPsi * sinPhi + cosTheta * sinPsi * cosPhi) * dy
			+ (sinPsi * sinTheta) * dz;

	ry = (-sinPsi * cosPhi - cosTheta * cosPsi * sinPhi) * dx
			+ (-sinPsi * sinPhi + cosTheta * cosPsi * cosPhi) * dy
			+ (cosPsi * sinTheta) * dz;

	rz = (sinTheta * sinPhi) * dx + (-sinTheta * cosPhi) * dy + (cosTheta) * dz;

	dx = rx;
	dy = ry;
	dz = rz;

	return *this;

} /* rotate */

//-*******************
// rotate(HepAxisAngle)
// rotate(HepEulerAngles)
//-*******************

//-***********************
// rotationOf(HepAxisAngle)
// rotationOf(HepEulerAngles)
// and coordinate axis rotations
//-***********************

Vector3 rotationOf(const Vector3 & vec, const Vector3 & axis, double delta) {
	Vector3 vv(vec);
	return vv.rotate(axis, delta);
}

Vector3 rotationOf(const Vector3 & vec, double phi, double theta, double psi) {
	Vector3 vv(vec);
	return vv.rotate(phi, theta, psi);
}

Vector3 rotationXOf(const Vector3 & vec, double delta) {
	Vector3 vv(vec);
	return vv.rotateX(delta);
}

Vector3 rotationYOf(const Vector3 & vec, double delta) {
	Vector3 vv(vec);
	return vv.rotateY(delta);
}

Vector3 rotationZOf(const Vector3 & vec, double delta) {
	Vector3 vv(vec);
	return vv.rotateZ(delta);
}



double Vector3::tolerance = Vector3::ToleranceTicks * 2.22045e-16;

} // namespace mpc

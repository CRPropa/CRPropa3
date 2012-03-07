// -*- C++ -*-
// $Id: ThreeVector.icc,v 1.2 2010/06/16 17:15:57 garren Exp $
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
// 
// This is the definitions of the inline member functions of the
// Vector3 class.
//

#include <cmath>

// ------------------
// Access to elements
// ------------------

// x, y, z
namespace mpc {
inline double & Vector3::operator[](int i) {
	return operator()(i);
}
inline double Vector3::operator[](int i) const {
	return operator()(i);
}

inline double Vector3::x() const {
	return dx;
}
inline double Vector3::y() const {
	return dy;
}
inline double Vector3::z() const {
	return dz;
}

inline double Vector3::getX() const {
	return dx;
}
inline double Vector3::getY() const {
	return dy;
}
inline double Vector3::getZ() const {
	return dz;
}

inline void Vector3::setX(double x) {
	dx = x;
}
inline void Vector3::setY(double y) {
	dy = y;
}
inline void Vector3::setZ(double z) {
	dz = z;
}

inline void Vector3::set(double x, double y, double z) {
	dx = x;
	dy = y;
	dz = z;
}

// --------------
// Global methods
// --------------

inline Vector3 Vector3::operator +(const Vector3 & b) const {
	return Vector3(x() + b.x(), y() + b.y(), z() + b.z());
}

inline Vector3 Vector3::operator -(const Vector3 & b) const {
	return Vector3(x() - b.x(), y() - b.y(), z() - b.z());
}

inline Vector3 Vector3::operator *( double a) const {
	return Vector3(a * x(), a * y(), a * z());
}

inline Vector3 operator *(double a, const Vector3 & p) {
	return Vector3(a * p.x(), a * p.y(), a * p.z());
}

inline double operator *(const Vector3 & a, const Vector3 & b) {
	return a.dot(b);
}

// --------------------------
// Set in various coordinates
// --------------------------

inline void Vector3::setRThetaPhi(double r, double theta, double phi) {
	setSpherical(r, theta, phi);
}

inline void Vector3::setREtaPhi(double r, double eta, double phi) {
	setSpherical(r, 2 * std::atan(std::exp(-eta)), phi);
}

inline void Vector3::setRhoPhiZ(double rho, double phi, double z) {
	setCylindrical(rho, phi, z);
}

// ------------
// Constructors
// ------------

inline Vector3::Vector3() :
	dx(0.), dy(0.), dz(0.) {
}
inline Vector3::Vector3(double x) :
	dx(x), dy(0.), dz(0.) {
}
inline Vector3::Vector3(double x, double y) :
	dx(x), dy(y), dz(0.) {
}
inline Vector3::Vector3(double x, double y, double z) :
	dx(x), dy(y), dz(z) {
}

inline Vector3::Vector3(const Vector3 & p) :
	dx(p.dx), dy(p.dy), dz(p.dz) {
}

inline Vector3::~Vector3() {
}

inline Vector3 & Vector3::operator =(const Vector3 & p) {
	dx = p.dx;
	dy = p.dy;
	dz = p.dz;
	return *this;
}

// ------------------
// Access to elements
// ------------------

// r, theta, phi

inline double Vector3::mag2() const {
	return dx * dx + dy * dy + dz * dz;
}
inline double Vector3::mag() const {
	return std::sqrt(mag2());
}
inline double Vector3::r() const {
	return mag();
}

inline double Vector3::theta() const {
	return dx == 0.0 && dy == 0.0 && dz == 0.0 ? 0.0 : std::atan2(perp(), dz);
}
inline double Vector3::phi() const {
	return dx == 0.0 && dy == 0.0 ? 0.0 : std::atan2(dy, dx);
}

inline double Vector3::getR() const {
	return mag();
}
inline double Vector3::getTheta() const {
	return theta();
}
inline double Vector3::getPhi() const {
	return phi();
}
inline double Vector3::angle() const {
	return theta();
}

inline double Vector3::cosTheta() const {
	double ptot = mag();
	return ptot == 0.0 ? 1.0 : dz / ptot;
}

inline double Vector3::cos2Theta() const {
	double ptot2 = mag2();
	return ptot2 == 0.0 ? 1.0 : dz * dz / ptot2;
}

inline void Vector3::setR(double r) {
	setMag(r);
}

inline void Vector3::setTheta(double th) {
	double ma = mag();
	double ph = phi();
	setX(ma * std::sin(th) * std::cos(ph));
	setY(ma * std::sin(th) * std::sin(ph));
	setZ(ma * std::cos(th));
}

inline void Vector3::setPhi(double ph) {
	double xy = perp();
	setX(xy * std::cos(ph));
	setY(xy * std::sin(ph));
}

// perp, eta, 

inline double Vector3::perp2() const {
	return dx * dx + dy * dy;
}
inline double Vector3::perp() const {
	return std::sqrt(perp2());
}
inline double Vector3::rho() const {
	return perp();
}
inline double Vector3::eta() const {
	return pseudoRapidity();
}

inline double Vector3::getRho() const {
	return perp();
}
inline double Vector3::getEta() const {
	return pseudoRapidity();
}

inline void Vector3::setPerp(double r) {
	double p = perp();
	if (p != 0.0) {
		dx *= r / p;
		dy *= r / p;
	}
}
inline void Vector3::setRho(double rho) {
	setPerp(rho);
}

// ----------
// Comparison
// ----------

inline bool Vector3::operator ==(const Vector3& v) const {
	return (v.x() == x() && v.y() == y() && v.z() == z()) ? true : false;
}

inline bool Vector3::operator !=(const Vector3& v) const {
	return (v.x() != x() || v.y() != y() || v.z() != z()) ? true : false;
}

inline double Vector3::getTolerance() {
	return tolerance;
}

// ----------
// Arithmetic
// ----------

inline Vector3& Vector3::operator +=(const Vector3 & p) {
	dx += p.x();
	dy += p.y();
	dz += p.z();
	return *this;
}

inline Vector3& Vector3::operator -=(const Vector3 & p) {
	dx -= p.x();
	dy -= p.y();
	dz -= p.z();
	return *this;
}

inline Vector3 Vector3::operator -() const {
	return Vector3(-dx, -dy, -dz);
}

inline Vector3& Vector3::operator *=(double a) {
	dx *= a;
	dy *= a;
	dz *= a;
	return *this;
}

// -------------------
// Combine two Vectors
// -------------------

inline double Vector3::diff2(const Vector3 & p) const {
	return (*this - p).mag2();
}

inline double Vector3::dot(const Vector3 & p) const {
	return dx * p.x() + dy * p.y() + dz * p.z();
}

inline Vector3 Vector3::cross(const Vector3 & p) const {
	return Vector3(dy * p.z() - p.y() * dz, dz * p.x() - p.z() * dx, dx
			* p.y() - p.x() * dy);
}

inline double Vector3::perp2(const Vector3 & p) const {
	double tot = p.mag2();
	double ss = dot(p);
	return tot > 0.0 ? mag2() - ss * ss / tot : mag2();
}

inline double Vector3::perp(const Vector3 & p) const {
	return std::sqrt(perp2(p));
}

inline Vector3 Vector3::perpPart() const {
	return Vector3(dx, dy, 0);
}
inline Vector3 Vector3::project() const {
	return Vector3(0, 0, dz);
}

inline Vector3 Vector3::perpPart(const Vector3 & v2) const {
	return (*this - project(v2));
}

inline double Vector3::angle(const Vector3 & q) const {
	return std::acos(cosTheta(q));
}

inline double Vector3::theta(const Vector3 & q) const {
	return angle(q);
}

inline double Vector3::azimAngle(const Vector3 & v2) const {
	return deltaPhi(v2);
}

// ----------
// Properties
// ----------

inline Vector3 Vector3::unit() const {
	double tot = mag2();
	Vector3 p(x(), y(), z());
	return tot > 0.0 ? p *= (1.0 / std::sqrt(tot)) : p;
}

inline Vector3 Vector3::orthogonal() const {
	double x = dx < 0.0 ? -dx : dx;
	double y = dy < 0.0 ? -dy : dy;
	double z = dz < 0.0 ? -dz : dz;
	if (x < y) {
		return x < z ? Vector3(0, dz, -dy) : Vector3(dy, -dx, 0);
	} else {
		return y < z ? Vector3(-dz, 0, dx) : Vector3(dy, -dx, 0);
	}
}

} // namespace mpc

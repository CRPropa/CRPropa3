// -*- C++ -*-
// CLASSDOC OFF
// $Id: ThreeVector.h,v 1.4 2010/06/16 17:15:57 garren Exp $
// ---------------------------------------------------------------------------
// CLASSDOC ON
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// Vector3 is a general 3-vector class defining vectors in three
// dimension using double components. Rotations of these vectors are
// performed by multiplying with an object of the HepRotation class.
//
// .SS See Also
// LorentzVector.h, Rotation.h, LorentzRotation.h 
//
// .SS Authors
// Leif Lonnblad and Anders Nilsson; Modified by Evgueni Tcherniaev;
// ZOOM additions by Mark Fischler
//

#ifndef HEP_THREEVECTOR_H
#define HEP_THREEVECTOR_H

#ifdef GNUPRAGMA
#pragma interface
#endif

#include <iostream>

namespace mpc {

/**
 @class Vector3
 @brief CLHEP 3-vector
 */
class Vector3 {

public:

	// Basic properties and operations on 3-vectors:

	enum {
		X = 0, Y = 1, Z = 2, NUM_COORDINATES = 3, SIZE = NUM_COORDINATES
	};
	// Safe indexing of the coordinates when using with matrices, arrays, etc.
	// (BaBar)

	Vector3();
	explicit Vector3(double x);
	Vector3(double x, double y);
	Vector3(double x, double y, double z);
	// The constructor.

	inline Vector3(const Vector3 &);
	// The copy constructor.

	inline ~Vector3();
	// The destructor.  Not virtual - inheritance from this class is dangerous.

	double operator ()(int) const;
	// Get components by index -- 0-based (Geant4)

	inline double operator [](int) const;
	// Get components by index -- 0-based (Geant4)

	double & operator ()(int);
	// Set components by index.  0-based.

	inline double & operator [](int);
	// Set components by index.  0-based.

	inline double x() const;
	inline double y() const;
	inline double z() const;
	// The components in cartesian coordinate system.  Same as getX() etc.

	inline void setX(double);
	inline void setY(double);
	inline void setZ(double);
	// Set the components in cartesian coordinate system.

	inline void set(double x, double y, double z);
	// Set all three components in cartesian coordinate system.

	inline double phi() const;
	// The azimuth angle.

	inline double theta() const;
	// The polar angle.

	inline double cosTheta() const;
	// Cosine of the polar angle.

	inline double cos2Theta() const;
	// Cosine squared of the polar angle - faster than cosTheta(). (ZOOM)

	inline double mag2() const;
	// The magnitude squared (r^2 in spherical coordinate system).

	inline double mag() const;
	// The magnitude (r in spherical coordinate system).

	inline void setPhi(double);
	// Set phi keeping mag and theta constant (BaBar).

	inline void setTheta(double);
	// Set theta keeping mag and phi constant (BaBar).

	void setMag(double);
	// Set magnitude keeping theta and phi constant (BaBar).

	inline double perp2() const;
	// The transverse component squared (rho^2 in cylindrical coordinate system).

	inline double perp() const;
	// The transverse component (rho in cylindrical coordinate system).

	inline void setPerp(double);
	// Set the transverse component keeping phi and z constant.

	void setCylTheta(double);
	// Set theta while keeping transvers component and phi fixed

	inline double perp2(const Vector3 &) const;
	// The transverse component w.r.t. given axis squared.

	inline double perp(const Vector3 &) const;
	// The transverse component w.r.t. given axis.

	inline Vector3 & operator =(const Vector3 &);
	// Assignment.

	inline bool operator ==(const Vector3 &) const;
	inline bool operator !=(const Vector3 &) const;
	// Comparisons (Geant4).

	bool isNear(const Vector3 &, double epsilon = tolerance) const;
	// Check for equality within RELATIVE tolerance (default 2.2E-14). (ZOOM)
	// |v1 - v2|**2 <= epsilon**2 * |v1.dot(v2)|

	double howNear(const Vector3 & v) const;
	// sqrt ( |v1-v2|**2 / v1.dot(v2) ) with a maximum of 1.
	// If v1.dot(v2) is negative, will return 1.

	double deltaR(const Vector3 & v) const;
	// sqrt( pseudorapity_difference**2 + deltaPhi **2 )

	inline Vector3 & operator +=(const Vector3 &);
	// Addition.

	inline Vector3 & operator -=(const Vector3 &);
	// Subtraction.

	inline Vector3 operator -() const;
	// Unary minus.

	inline Vector3 & operator *=(double);
	// Scaling with real numbers.

	Vector3 & operator /=(double);
	// Division by (non-zero) real number.

	inline Vector3 unit() const;
	// Vector parallel to this, but of length 1.

	inline Vector3 orthogonal() const;
	// Vector orthogonal to this (Geant4).

	inline double dot(const Vector3 &) const;
	// double product.

	inline Vector3 cross(const Vector3 &) const;
	// Cross product.

	double angle(const Vector3 &) const;
	// The angle w.r.t. another 3-vector.

	double pseudoRapidity() const;
	// Returns the pseudo-rapidity, i.e. -ln(tan(theta/2))

	void setEta(double p);
	// Set pseudo-rapidity, keeping magnitude and phi fixed.  (ZOOM)

	void setCylEta(double p);
	// Set pseudo-rapidity, keeping transverse component and phi fixed.  (ZOOM)

	Vector3 & rotateX(double);
	// Rotates the Vector3 around the x-axis.

	Vector3 & rotateY(double);
	// Rotates the Vector3 around the y-axis.

	Vector3 & rotateZ(double);
	// Rotates the Vector3 around the z-axis.

	Vector3 & rotateUz(const Vector3&);
	// Rotates reference frame from Uz to newUz (unit vector) (Geant4).

	//Vector3 & rotate(double, const Vector3 &);
	// Rotates around the axis specified by another Vector3.
	// (Uses methods of HepRotation, forcing linking in of Rotation.cc.)

	// Transformation with a Rotation matrix.

	// = = = = = = = = = = = = = = = = = = = = = = = =
	//
	// Esoteric properties and operations on 3-vectors:
	//
	// 1 - Set vectors in various coordinate systems
	// 2 - Synonyms for accessing coordinates and properties
	// 3 - Comparisions (dictionary, near-ness, and geometric)
	// 4 - Intrinsic properties
	// 5 - Properties releative to z axis and arbitrary directions
	// 6 - Polar and azimuthal angle decomposition and deltaPhi
	// 7 - Rotations
	//
	// = = = = = = = = = = = = = = = = = = = = = = = =

	// 1 - Set vectors in various coordinate systems

	inline void setRThetaPhi(double r, double theta, double phi);
	// Set in spherical coordinates:  Angles are measured in RADIANS

	inline void setREtaPhi(double r, double eta, double phi);
	// Set in spherical coordinates, but specify peudorapidiy to determine theta.

	inline void setRhoPhiZ(double rho, double phi, double z);
	// Set in cylindrical coordinates:  Phi angle is measured in RADIANS

	void setRhoPhiTheta(double rho, double phi, double theta);
	// Set in cylindrical coordinates, but specify theta to determine z.

	void setRhoPhiEta(double rho, double phi, double eta);
	// Set in cylindrical coordinates, but specify pseudorapidity to determine z.

	// 2 - Synonyms for accessing coordinates and properties

	inline double getX() const;
	inline double getY() const;
	inline double getZ() const;
	// x(), y(), and z()

	inline double getR() const;
	inline double getTheta() const;
	inline double getPhi() const;
	// mag(), theta(), and phi()

	inline double r() const;
	// mag()

	inline double rho() const;
	inline double getRho() const;
	// perp()

	double eta() const;
	double getEta() const;
	// pseudoRapidity()

	inline void setR(double s);
	// setMag()

	inline void setRho(double s);
	// setPerp()

	// 3 - Comparisions (dictionary, near-ness, and geometric)

	int compare(const Vector3 & v) const;
	bool operator >(const Vector3 & v) const;
	bool operator <(const Vector3 & v) const;
	bool operator>=(const Vector3 & v) const;
	bool operator<=(const Vector3 & v) const;
	// dictionary ordering according to z, then y, then x component

	inline double diff2(const Vector3 & v) const;
	// |v1-v2|**2

	static double setTolerance(double tol);
	static inline double getTolerance();
	// Set the tolerance used in isNear() for Vector3s

	bool isParallel(const Vector3 & v, double epsilon = tolerance) const;
	// Are the vectors parallel, within the given tolerance?

	bool isOrthogonal(const Vector3 & v, double epsilon = tolerance) const;
	// Are the vectors orthogonal, within the given tolerance?

	double howParallel(const Vector3 & v) const;
	// | v1.cross(v2) / v1.dot(v2) |, to a maximum of 1.

	double howOrthogonal(const Vector3 & v) const;
	// | v1.dot(v2) / v1.cross(v2) |, to a maximum of 1.

	enum {
		ToleranceTicks = 100
	};

	// 4 - Intrinsic properties

	double beta() const;
	// relativistic beta (considering v as a velocity vector with c=1)
	// Same as mag() but will object if >= 1

	double gamma() const;
	// relativistic gamma (considering v as a velocity vector with c=1)

	double coLinearRapidity() const;
	// inverse tanh (beta)

	// 5 - Properties relative to Z axis and to an arbitrary direction

	// Note that the non-esoteric CLHEP provides
	// theta(), cosTheta(), cos2Theta, and angle(const Vector3&)

	inline double angle() const;
	// angle against the Z axis -- synonym for theta()

	inline double theta(const Vector3 & v2) const;
	// synonym for angle(v2)

	double cosTheta(const Vector3 & v2) const;
	double cos2Theta(const Vector3 & v2) const;
	// cos and cos^2 of the angle between two vectors

	inline Vector3 project() const;
	Vector3 project(const Vector3 & v2) const;
	// projection of a vector along a direction.

	inline Vector3 perpPart() const;
	inline Vector3 perpPart(const Vector3 & v2) const;
	// vector minus its projection along a direction.

	double rapidity() const;
	// inverse tanh(v.z())

	double rapidity(const Vector3 & v2) const;
	// rapidity with respect to specified direction:
	// inverse tanh (v.dot(u)) where u is a unit in the direction of v2

	double eta(const Vector3 & v2) const;
	// - ln tan of the angle beween the vector and the ref direction.

	// 6 - Polar and azimuthal angle decomposition and deltaPhi

	// Decomposition of an angle within reference defined by a direction:

	double polarAngle(const Vector3 & v2) const;
	// The reference direction is Z: the polarAngle is abs(v.theta()-v2.theta()).

	double deltaPhi(const Vector3 & v2) const;
	// v.phi()-v2.phi(), brought into the range (-PI,PI]

	double azimAngle(const Vector3 & v2) const;
	// The reference direction is Z: the azimAngle is the same as deltaPhi

	double polarAngle(const Vector3 & v2, const Vector3 & ref) const;
	// For arbitrary reference direction,
	// 	polarAngle is abs(v.angle(ref) - v2.angle(ref)).

	double azimAngle(const Vector3 & v2, const Vector3 & ref) const;
	// To compute azimangle, project v and v2 into the plane normal to
	// the reference direction.  Then in that plane take the angle going
	// clockwise around the direction from projection of v to that of v2.

	// 7 - Rotations

	// These mehtods **DO NOT** use anything in the HepRotation class.
	// Thus, use of v.rotate(axis,delta) does not force linking in Rotation.cc.

	Vector3 & rotate(const Vector3 & axis, double delta);
	// Synonym for rotate (delta, axis)

	Vector3 & rotate(double phi, double theta, double psi);
	// Rotate via Euler Angles. Our Euler Angles conventions are
	// those of Goldstein Classical Mechanics page 107.

	inline Vector3 operator +(const Vector3 &) const;
	// Addition of 3-vectors.

	inline Vector3 operator -(const Vector3 &) const;
	Vector3 operator /(double a) const;
	inline Vector3 operator *( double a)const;

protected:
	void setSpherical(double r, double theta, double phi);
	void setCylindrical(double r, double phi, double z);
	double negativeInfinity() const;

protected:

	double dx;
	double dy;
	double dz;
	// The components.

	static double tolerance;
	// default tolerance criterion for isNear() to return true.
};
// Vector3

// Global Methods
//
//Vector3 rotationXOf(const Vector3 & vec, double delta);
//Vector3 rotationYOf(const Vector3 & vec, double delta);
//Vector3 rotationZOf(const Vector3 & vec, double delta);

//Vector3 rotationOf(const Vector3 & vec, const Vector3 & axis, double delta);
//
//Vector3 rotationOf(const Vector3 & vec, double phi, double theta, double psi);

// Return a new vector based on a rotation of the supplied vector

std::ostream & operator <<(std::ostream &, const Vector3 &);
// Output to a stream.

std::istream & operator >>(std::istream &, Vector3 &);
// Input from a stream.

//extern const Vector3 HepXHat, HepYHat, HepZHat;
//
//typedef Vector3 HepThreeVectorD;
//typedef Vector3 HepThreeVectorF;
//
// Division of 3-vectors by non-zero real number


// Subtraction of 3-vectors.

inline double operator *(const Vector3 &, const Vector3 &);
// double product of 3-vectors.

inline Vector3 operator *(double a, const Vector3 &);
//// Scaling of 3-vectors with a real number

}// namespace mpc

#include "mpc/Vector3Inline.h"

#endif /* HEP_THREEVECTOR_H */

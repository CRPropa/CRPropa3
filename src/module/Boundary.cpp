#include "crpropa/module/Boundary.h"
#include "crpropa/Units.h"

#include <sstream>

namespace crpropa {

PeriodicBox::PeriodicBox() :
		origin(Vector3d(0, 0, 0)), size(Vector3d(0, 0, 0)) {
}

PeriodicBox::PeriodicBox(Vector3d o, Vector3d s) :
		origin(o), size(s) {
}

void PeriodicBox::process(Candidate *c) const {
	Vector3d pos = c->current.getPosition();
	Vector3d n = ((pos - origin) / size).floor();

	if ((n.x == 0) and (n.y == 0) and (n.z == 0))
		return; // do nothing if candidate is inside the box

	c->current.setPosition(pos - n * size);
	c->previous.setPosition(c->previous.getPosition() - n * size);
	c->source.setPosition(c->source.getPosition() - n * size);
	c->created.setPosition(c->created.getPosition() - n * size);
}

void PeriodicBox::setOrigin(Vector3d o) {
	origin = o;
}
void PeriodicBox::setSize(Vector3d s) {
	size = s;
}

std::string PeriodicBox::getDescription() const {
	std::stringstream s;
	s << "Periodic box: origin " << origin / Mpc << " Mpc, ";
	s << "size " << size / Mpc << " Mpc";
	return s.str();
}


ReflectiveShell::ReflectiveShell(Vector3d center, double r) :
		center(center), radius(r) {
		}

double ReflectiveShell::distance(const Vector3d &point) const {
	Vector3d dR = point - center;
	return dR.getR() - radius;
}

Vector3d ReflectiveShell::normal(const Vector3d& point) const {
	Vector3d d = point - center;
	return d.getUnitVector();
}

void ReflectiveShell::process(Candidate *c) const {
	double currentDistance = distance(c->current.getPosition());
	double previousDistance = distance(c->previous.getPosition());
	// check if cosmic ray crossed boundary in last step
	if (currentDistance * previousDistance < 0){
		Vector3d currentDirection = c->current.getDirection();
		Vector3d previousPosition = c->previous.getPosition();
		// get point where trajectory intersects the shell boundary
		double p_half = previousPosition.dot(currentDirection) / currentDirection.getR2();
		double q = (previousPosition.getR2() - radius * radius) / currentDirection.getR2();
		double k1 = - p_half + sqrt(p_half * p_half - q);
		Vector3d intersectPoint = previousPosition + k1 * currentDirection;
		// flip component of velocity normal to surface
		Vector3d surfaceNormal = normal(intersectPoint);
		Vector3d newDirection = currentDirection - 2 * currentDirection.dot(surfaceNormal) * surfaceNormal;
		c->current.setDirection(newDirection.getUnitVector());
		// also update position
		c->current.setPosition(intersectPoint + newDirection * (c->getCurrentStep() - (k1 * currentDirection).getR()) / newDirection.getR());
	}
}

void ReflectiveShell::setCenter(Vector3d c) {
	center = c;
}
void ReflectiveShell::setRadius(double r) {
	radius = r;
}

std::string ReflectiveShell::getDescription() const {
	std::stringstream s;
	s << "Reflective sphere: center: " << center / Mpc << " Mpc, ";
	s << "radius: " << radius / Mpc << " Mpc";
	return s.str();
};


ReflectiveBox::ReflectiveBox() :
		origin(Vector3d(0, 0, 0)), size(Vector3d(0, 0, 0)) {
}

ReflectiveBox::ReflectiveBox(Vector3d o, Vector3d s) :
		origin(o), size(s) {
}

void ReflectiveBox::process(Candidate *c) const {
	Vector3d cur = (c->current.getPosition() - origin) / size; // current position in cell units
	Vector3d n = cur.floor();

	if ((n.x == 0) and (n.y == 0) and (n.z == 0))
		return; // do nothing if candidate is inside the box

	// flip direction
	Vector3d nReflect(pow(-1, n.x), pow(-1, n.y), pow(-1, n.z));
	c->current.setDirection(c->current.getDirection() * nReflect);
	c->previous.setDirection(c->previous.getDirection() * nReflect);
	c->created.setDirection(c->created.getDirection() * nReflect);
	c->source.setDirection(c->source.getDirection() * nReflect);

	Vector3d src = (c->source.getPosition() - origin) / size; // initial position in cell units
	Vector3d cre = (c->created.getPosition() - origin) / size; // initial position in cell units
	Vector3d prv = (c->previous.getPosition() - origin) / size; // previous position in cell units

	// repeatedly translate until the current position is inside the cell
	while ((cur.x < 0) or (cur.x > 1)) {
		double t = 2 * (cur.x > 1);
		src.x = t - src.x;
		cre.x = t - cre.x;
		prv.x = t - prv.x;
		cur.x = t - cur.x;
	}
	while ((cur.y < 0) or (cur.y > 1)) {
		double t = 2 * (cur.y > 1);
		src.y = t - src.y;
		cre.y = t - cre.y;
		prv.y = t - prv.y;
		cur.y = t - cur.y;
	}
	while ((cur.z < 0) or (cur.z > 1)) {
		double t = 2 * (cur.z > 1);
		src.z = t - src.z;
		cre.z = t - cre.z;
		prv.z = t - prv.z;
		cur.z = t - cur.z;
	}

	c->current.setPosition(cur * size + origin);
	c->source.setPosition(src * size + origin);
	c->created.setPosition(cre * size + origin);
	c->previous.setPosition(prv * size + origin);
}

void ReflectiveBox::setOrigin(Vector3d o) {
	origin = o;
}
void ReflectiveBox::setSize(Vector3d s) {
	size = s;
}

std::string ReflectiveBox::getDescription() const {
	std::stringstream s;
	s << "Reflective box: origin " << origin / Mpc << " Mpc, ";
	s << "size " << size / Mpc << " Mpc";
	return s.str();
}

CubicBoundary::CubicBoundary() :
		origin(Vector3d(0, 0, 0)), size(0), limitStep(true), margin(0.1 * kpc) {
}

CubicBoundary::CubicBoundary(Vector3d o, double s) :
		origin(o), size(s), limitStep(true), margin(0.1 * kpc) {
}

void CubicBoundary::process(Candidate *c) const {
	Vector3d r = c->current.getPosition() - origin;
	double lo = r.min();
	double hi = r.max();
	if ((lo <= 0) or (hi >= size)) {
		reject(c);
	}
	if (limitStep) {
		c->limitNextStep(lo + margin);
		c->limitNextStep(size - hi + margin);
	}
}

void CubicBoundary::setOrigin(Vector3d o) {
	origin = o;
}
void CubicBoundary::setSize(double s) {
	size = s;
}
void CubicBoundary::setMargin(double m) {
	margin = m;
}
void CubicBoundary::setLimitStep(bool b) {
	limitStep = b;
}

std::string CubicBoundary::getDescription() const {
	std::stringstream s;
	s << "Cubic Boundary: origin " << origin / Mpc << " Mpc, ";
	s << "size " << size / Mpc << " Mpc, ";
	s << "Flag: '" << rejectFlagKey << "' -> '" << rejectFlagValue << "', ";
	s << "MakeInactive: " << (makeRejectedInactive ? "yes" : "no");
	if (rejectAction.valid())
		s << ", Action: " << rejectAction->getDescription();
	return s.str();
}

SphericalBoundary::SphericalBoundary() :
		center(Vector3d(0, 0, 0)), radius(0), limitStep(true), margin(0.1 * kpc) {
}

SphericalBoundary::SphericalBoundary(Vector3d c, double r) :
		center(c), radius(r), limitStep(true), margin(0.1 * kpc) {
}

void SphericalBoundary::process(Candidate *c) const {
	double d = (c->current.getPosition() - center).getR();
	if (d >= radius) {
		reject(c);
	}
	if (limitStep)
		c->limitNextStep(radius - d + margin);
}

void SphericalBoundary::setCenter(Vector3d c) {
	center = c;
}
void SphericalBoundary::setRadius(double r) {
	radius = r;
}
void SphericalBoundary::setMargin(double m) {
	margin = m;
}
void SphericalBoundary::setLimitStep(bool b) {
	limitStep = b;
}

std::string SphericalBoundary::getDescription() const {
	std::stringstream s;
	s << "Spherical Boundary: radius " << radius / Mpc << " Mpc, ";
	s << "around " << center / Mpc << " Mpc, ";
	s << "Flag: '" << rejectFlagKey << "' -> '" << rejectFlagValue << "', ";
	s << "MakeInactive: " << (makeRejectedInactive ? "yes" : "no");
	if (rejectAction.valid())
		s << ", Action: " << rejectAction->getDescription();
	return s.str();
}

EllipsoidalBoundary::EllipsoidalBoundary() :
		focalPoint1(Vector3d(0, 0, 0)), focalPoint2(Vector3d(0, 0, 0)),
		majorAxis(0), limitStep(true), margin(0.1 * kpc) {
}

EllipsoidalBoundary::EllipsoidalBoundary(Vector3d f1, Vector3d f2, double a) :
		focalPoint1(f1), focalPoint2(f2), majorAxis(a), limitStep(true),
		margin(0.1 * kpc) {
}

void EllipsoidalBoundary::process(Candidate *c) const {
	Vector3d pos = c->current.getPosition();
	double d = pos.getDistanceTo(focalPoint1) + pos.getDistanceTo(focalPoint2);
	if (d >= majorAxis) {
		reject(c);
	}
	if (limitStep)
		c->limitNextStep(majorAxis - d + margin);
}

void EllipsoidalBoundary::setFocalPoints(Vector3d f1, Vector3d f2) {
	focalPoint1 = f1;
	focalPoint2 = f2;
}
void EllipsoidalBoundary::setMajorAxis(double a) {
	majorAxis = a;
}
void EllipsoidalBoundary::setMargin(double m) {
	margin = m;
}
void EllipsoidalBoundary::setLimitStep(bool b) {
	limitStep = b;
}

std::string EllipsoidalBoundary::getDescription() const {
	std::stringstream s;
	s << "Ellipsoidal Boundary: F1 = " << focalPoint1 / Mpc << " Mpc, ";
	s << "F2 = " << focalPoint2 / Mpc << " Mpc, ";
	s << "major axis " << majorAxis / Mpc << " Mpc; ";
	s << "Flag: '" << rejectFlagKey << "' -> '" << rejectFlagValue << "', ";
	s << "MakeInactive: " << (makeRejectedInactive ? "yes" : "no");
	if (rejectAction.valid())
		s << ", Action: " << rejectAction->getDescription();
	return s.str();
}

CylindricalBoundary::CylindricalBoundary() :
  origin(Vector3d(0,0,0)), height(0), radius(0), limitStep(false), margin(0) {
}

CylindricalBoundary::CylindricalBoundary(Vector3d o, double h, double r) :
  origin(o), height(h), radius(r), limitStep(false) , margin(0){
}

void CylindricalBoundary::process(Candidate *c) const {
	Vector3d d = c->current.getPosition() - origin;
	double R2 = pow(d.x, 2.)+pow(d.y, 2.);
	double Z = fabs(d.z);
	if ( R2 < pow(radius, 2.) and Z < height/2.) {
		if(limitStep) {
			c->limitNextStep(std::min(radius - pow(R2, 0.5), height/2. - Z) + margin);	
		}
	  return;
	}
	reject(c);
}

void CylindricalBoundary::setOrigin(Vector3d o) {
	origin = o;
}
void CylindricalBoundary::setHeight(double h) {
	height = h;
}
void CylindricalBoundary::setRadius(double r) {
	radius = r;
}
void CylindricalBoundary::setLimitStep(bool b) {
	limitStep = b;
}
void CylindricalBoundary::setMargin(double m) {
	margin = m;
}


std::string CylindricalBoundary::getDescription() const {
	std::stringstream s;
	s << "Cylindrical Boundary: origin = " << origin / kpc << " kpc, ";
	s << "radius = " << radius / kpc << " kpc, ";
	s << "height" << height / kpc << " kpc; ";
	s << "Flag: '" << rejectFlagKey << "' -> '" << rejectFlagValue << "', ";
	s << "MakeInactive: " << (makeRejectedInactive ? "yes" : "no");
	if (rejectAction.valid())
		s << ", Action: " << rejectAction->getDescription();
	return s.str();
}

} // namespace crpropa

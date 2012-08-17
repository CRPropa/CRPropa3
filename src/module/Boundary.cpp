#include "mpc/module/Boundary.h"

#include <sstream>

namespace mpc {

PeriodicBox::PeriodicBox(Vector3d o, Vector3d s) :
		origin(o), size(s) {
	updateDescription();
}

void PeriodicBox::process(Candidate *c) const {
	Vector3d pos = c->current.getPosition();
	Vector3d n = ((pos - origin) / size).floor();

	if ((n.x == 0) and (n.y == 0) and (n.z == 0))
		return; // do nothing if candidate is inside the box

	c->initial.setPosition(c->initial.getPosition() - n * size);
	c->current.setPosition(pos - n * size);
}

void PeriodicBox::updateDescription() {
	std::stringstream s;
	s << "Periodic box: origin " << origin << ", size " << size;
	setDescription(s.str());
}

ReflectiveBox::ReflectiveBox(Vector3d o, Vector3d s) :
		origin(o), size(s) {
	updateDescription();
}

void ReflectiveBox::process(Candidate *c) const {
	Vector3d cur = (c->current.getPosition() - origin) / size; // current position in unit cells
	Vector3d n = cur.floor();

	if ((n.x == 0) and (n.y == 0) and (n.z == 0))
		return; // do nothing if candidate is inside the box

	// flip direction
	Vector3d nReflect(pow(-1, n.x), pow(-1, n.y), pow(-1, n.z));
	c->current.setDirection(c->current.getDirection() * nReflect);
	c->initial.setDirection(c->initial.getDirection() * nReflect);

	Vector3d ini = (c->initial.getPosition() - origin) / size; // initial position in unit cells

	// repeatedly translate until the current position is inside the cell
	while ((cur.x < 0) or (cur.x > 1)) {
		ini.x = 2 * (cur.x > 1) - ini.x;
		cur.x = 2 * (cur.x > 1) - cur.x;
	}
	while ((cur.y < 0) or (cur.y > 1)) {
		ini.y = 2 * (cur.y > 1) - ini.y;
		cur.y = 2 * (cur.y > 1) - cur.y;
	}
	while ((cur.z < 0) or (cur.z > 1)) {
		ini.z = 2 * (cur.z > 1) - ini.z;
		cur.z = 2 * (cur.z > 1) - cur.z;
	}

	c->current.setPosition(cur * size + origin);
	c->initial.setPosition(ini * size + origin);
}

void ReflectiveBox::updateDescription() {
	std::stringstream s;
	s << "Reflective box: origin " << origin << ", size " << size;
	setDescription(s.str());
}

CubicBoundary::CubicBoundary(Vector3d o, double s, std::string f, std::string v) :
		origin(o), size(s), margin(0), flag(f), flagValue(v), limitStep(false) {
	updateDescription();
}

void CubicBoundary::setLimitStep(bool limitStep, double margin) {
	this->limitStep = limitStep;
	this->margin = margin;
	updateDescription();
}

void CubicBoundary::process(Candidate *c) const {
	Vector3d r = c->current.getPosition() - origin;
	double lo = r.min();
	double hi = r.max();
	if ((lo <= 0) or (hi >= size)) {
		c->setActive(false);
		c->setProperty(flag, flagValue);
	}
	if (limitStep) {
		c->limitNextStep(lo + margin);
		c->limitNextStep(size - hi + margin);
	}
}

void CubicBoundary::updateDescription() {
	std::stringstream s;
	s << "Cubic Boundary: origin " << origin << ", size " << size;
	s << " Flag: " << flag << " -> " << flagValue;
	setDescription(s.str());
}

SphericalBoundary::SphericalBoundary(Vector3d c, double r, std::string f,
		std::string v) :
		center(c), radius(r), flag(f), flagValue(v), limitStep(false), margin(0) {
	updateDescription();
}

void SphericalBoundary::setLimitStep(bool limitStep, double margin) {
	this->limitStep = limitStep;
	this->margin = margin;
	updateDescription();
}

void SphericalBoundary::process(Candidate *c) const {
	double d = (c->current.getPosition() - center).getMag();
	if (d >= radius) {
		c->setActive(false);
		c->setProperty(flag, flagValue);
	}
	if (limitStep)
		c->limitNextStep(radius - d + margin);
}

void SphericalBoundary::updateDescription() {
	std::stringstream s;
	s << "Spherical Boundary: radius " << radius << " around " << center;
	s << " Flag: " << flag << " -> " << flagValue;
	setDescription(s.str());
}

EllipsoidalBoundary::EllipsoidalBoundary(Vector3d f1,
		Vector3d f2, double a, std::string f,
		std::string v) : focalPoint1(f1), focalPoint2(f2), majorAxis(a), flag(f), flagValue(v), limitStep(false), margin(0) {
	updateDescription();
}

void EllipsoidalBoundary::setLimitStep(bool limitStep, double margin) {
	this->limitStep = limitStep;
	this->margin = margin;
	updateDescription();
}

void EllipsoidalBoundary::process(Candidate *c) const {
	Vector3d pos = c->current.getPosition();
	double d = pos.getDistanceTo(focalPoint1) + pos.getDistanceTo(focalPoint2);
	if (d >= majorAxis) {
		c->setActive(false);
		c->setProperty(flag, flagValue);
	}
	if (limitStep)
		c->limitNextStep(majorAxis - d + margin);
}

void EllipsoidalBoundary::updateDescription() {
	std::stringstream s;
	s << "Ellipsoidal Boundary: F1 = " << focalPoint1 / Mpc << ", F2 = "
			<< focalPoint2 / Mpc << ", major axis = " << majorAxis / Mpc;
	s << " Flag: " << flag << " -> " << flagValue;
	setDescription(s.str());
}

} // namespace mpc

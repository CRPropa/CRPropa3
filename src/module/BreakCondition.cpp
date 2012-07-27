#include "mpc/module/BreakCondition.h"

namespace mpc {

MaximumTrajectoryLength::MaximumTrajectoryLength(double maxLength) {
	this->maxLength = maxLength;
	updateDescription();
}

void MaximumTrajectoryLength::setMaximumTrajectoryLength(double length) {
	maxLength = length;
	updateDescription();
}

double MaximumTrajectoryLength::getMaximumTrajectoryLength() const {
	return maxLength;
}

void MaximumTrajectoryLength::process(Candidate *c) const {
	double l = c->getTrajectoryLength();
	if (l >= maxLength) {
		c->setActive(false);
		c->setProperty("Deactivated", getDescription());
	} else {
		c->limitNextStep(maxLength - l);
	}
}

void MaximumTrajectoryLength::updateDescription() {
	std::stringstream s;
	s << "Maximum trajectory length: " << maxLength / Mpc << " Mpc";
	setDescription(s.str());
}

MinimumEnergy::MinimumEnergy(double minEnergy) {
	this->minEnergy = minEnergy;
	updateDescription();
}

void MinimumEnergy::setMinimumEnergy(double energy) {
	minEnergy = energy;
	updateDescription();
}

double MinimumEnergy::getMinimumEnergy() const {
	return minEnergy;
}

void MinimumEnergy::process(Candidate *c) const {
	if (c->current.getEnergy() <= minEnergy) {
		c->setActive(false);
		c->setProperty("Deactivated", getDescription());
	}
}

void MinimumEnergy::updateDescription() {
	std::stringstream s;
	s << "Minimum energy: " << minEnergy / EeV << " EeV";
	setDescription(s.str());
}

SmallObserverSphere::SmallObserverSphere(Vector3d center, double radius,
		std::string flag, std::string flagValue) {
	this->center = center;
	this->radius = radius;
	this->flag = flag;
	this->flagValue = flagValue;
	this->makeInactive = true;
	updateDescription();
}

void SmallObserverSphere::setMakeInactive(bool makeInactive) {
	this->makeInactive = makeInactive;
	updateDescription();
}

void SmallObserverSphere::process(Candidate *c) const {
	double d = (c->current.getPosition() - center).getMag();
	if (d <= radius * 1.01) {
		c->setProperty(flag, flagValue);
		if (makeInactive) {
			c->setActive(false);
			c->setProperty("Deactivated", getDescription());
		}
	}
	c->limitNextStep((d - radius));
}

void SmallObserverSphere::updateDescription() {
	std::stringstream s;
	s << "Small observer sphere: " << radius / Mpc << " Mpc radius around "
			<< center / Mpc;
	s << " Mpc, Flag: '" << flag << "' -> '" << flagValue << "'";
	setDescription(s.str());
}

PeriodicBox::PeriodicBox(Vector3d origin, Vector3d size) {
	this->origin = origin;
	this->size = size;
	updateDescription();
}

void PeriodicBox::process(Candidate *c) const {
	Vector3d relPos = c->current.getPosition() - origin;
	Vector3d n = (relPos / size).floor(); // integers for translation
	if ((n.x != 0) or (n.y != 0) or (n.z != 0)) {
		c->initial.setPosition(c->initial.getPosition() - n * size);
		c->current.setPosition(c->current.getPosition() - n * size);
	}
}

void PeriodicBox::updateDescription() {
	std::stringstream s;
	s << "Periodic box: origin " << origin << ", size " << size;
	setDescription(s.str());
}

ReflectiveBox::ReflectiveBox(Vector3d origin, Vector3d size) {
	this->origin = origin;
	this->size = size;
	updateDescription();
}

void ReflectiveBox::process(Candidate *c) const {
	Vector3d relPos = c->current.getPosition() - origin;
	Vector3d n = (relPos / size).floor(); // integers for translation
	std::cout << relPos << std::endl;
	std::cout << n << std::endl;
	if ((n.x != 0) or (n.y != 0) or (n.z != 0)) {
		// flip direction
		Vector3d dir = c->current.getDirection();
		dir.x *= pow(-1, n.x);
		dir.y *= pow(-1, n.y);
		dir.z *= pow(-1, n.z);
		c->current.setDirection(dir);
		// translate into the box
		Vector3d t = (relPos % size) * n;
		c->initial.setPosition(c->initial.getPosition() + t);
		c->current.setPosition(c->current.getPosition() + t);
	}
}

void ReflectiveBox::updateDescription() {
	std::stringstream s;
	s << "Reflective box: origin " << origin << ", size " << size;
	setDescription(s.str());
}

CubicBoundary::CubicBoundary(Vector3d origin, double size, std::string flag,
		std::string flagValue) {
	this->origin = origin;
	this->size = size;
	this->flag = flag;
	this->flagValue = flagValue;
	this->makeInactive = true;
	this->limitStep = false;
	this->margin = 0;
	updateDescription();
}

void CubicBoundary::setMakeInactive(bool makeInactive) {
	this->makeInactive = makeInactive;
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
		c->setProperty(flag, flagValue);
		if (makeInactive) {
			c->setActive(false);
			c->setProperty("Deactivated", getDescription());
		}
	} else if (limitStep) {
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

SphericalBoundary::SphericalBoundary(Vector3d center, double radius,
		std::string flag, std::string flagValue) {
	this->center = center;
	this->radius = radius;
	this->flag = flag;
	this->flagValue = flagValue;
	this->makeInactive = true;
	this->limitStep = false;
	this->margin = 0;
	updateDescription();
}

void SphericalBoundary::setMakeInactive(bool makeInactive) {
	this->makeInactive = makeInactive;
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
		c->setProperty(flag, flagValue);
		if (makeInactive) {
			c->setActive(false);
			c->setProperty("Deactivated", getDescription());
		}
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

EllipsoidalBoundary::EllipsoidalBoundary(Vector3d focalPoint1,
		Vector3d focalPoint2, double majorAxis, std::string flag,
		std::string flagValue) {
	this->focalPoint1 = focalPoint1;
	this->focalPoint2 = focalPoint2;
	this->majorAxis = majorAxis;
	this->flag = flag;
	this->flagValue = flagValue;
	this->makeInactive = true;
	this->limitStep = false;
	this->margin = 0;
	updateDescription();
}

void EllipsoidalBoundary::setMakeInactive(bool makeInactive) {
	this->makeInactive = makeInactive;
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
		c->setProperty(flag, flagValue);
		if (makeInactive) {
			c->setActive(false);
			c->setProperty("Deactivated", getDescription());
		}
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

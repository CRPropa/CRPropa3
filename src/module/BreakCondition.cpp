#include "mpc/module/BreakCondition.h"

namespace mpc {

MaximumTrajectoryLength::MaximumTrajectoryLength(double length) :
		maxLength(length) {
	std::stringstream s;
	s << "Maximum trajectory length: " << maxLength / Mpc << " Mpc";
	setDescription(s.str());
}

void MaximumTrajectoryLength::process(Candidate *candidate) const {
	double l = candidate->getTrajectoryLength();
	if (l >= maxLength) {
		candidate->setActive(false);
		candidate->setProperty("Deactivated", getDescription());
	} else
		candidate->limitNextStep(maxLength - l);
}

MinimumEnergy::MinimumEnergy(double minEnergy) :
		minEnergy(minEnergy) {
	std::stringstream s;
	s << "Minimum energy: " << minEnergy / EeV << " EeV";
	setDescription(s.str());
}

void MinimumEnergy::process(Candidate *candidate) const {
	if (candidate->current.getEnergy() <= minEnergy) {
		candidate->setActive(false);
		candidate->setProperty("Deactivated", getDescription());
	}
}

SmallObserverSphere::SmallObserverSphere(Vector3 center, double radius,
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

void SmallObserverSphere::process(Candidate *candidate) const {
	double d = (candidate->current.getPosition() - center).mag();
	if (d <= radius * 1.01) {
		candidate->setProperty("Detected", "");
		if (makeInactive) {
			candidate->setActive(false);
			candidate->setProperty("Deactivated", getDescription());
		}
	}
	candidate->limitNextStep((d - radius));
}

void SmallObserverSphere::updateDescription() {
	std::stringstream s;
	s << "Small observer sphere: " << radius << " radius around " << center;
	s << " Flag: " << flag << " -> " << flagValue;
	setDescription(s.str());
}

CubicBoundary::CubicBoundary(Vector3 origin, double size, std::string flag,
		std::string flagValue) {
	this->origin = origin;
	this->size = size;
	this->flag = flag;
	this->flagValue = flagValue;
	this->makeInactive = false;
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

void CubicBoundary::process(Candidate *candidate) const {
	Vector3 relPos = candidate->current.getPosition() - origin;
	double lo = std::min(relPos.x(), std::min(relPos.y(), relPos.z()));
	double hi = std::max(relPos.x(), std::max(relPos.y(), relPos.z()));
	if ((lo <= 0.) or (hi >= size)) {
		candidate->setProperty(flag, flagValue);
		if (makeInactive) {
			candidate->setActive(false);
			candidate->setProperty("Deactivated", getDescription());
		}
	}
	if (limitStep) {
		candidate->limitNextStep(lo + margin);
		candidate->limitNextStep(size - hi + margin);
	}
}

void CubicBoundary::updateDescription() {
	std::stringstream s;
	s << "Cubic Boundary: origin " << origin << ", size " << size;
	s << " Flag: " << flag << " -> " << flagValue;
	setDescription(s.str());
}

SphericalBoundary::SphericalBoundary(Vector3 center, double radius,
		std::string flag, std::string flagValue) {
	this->center = center;
	this->radius = radius;
	this->flag = flag;
	this->flagValue = flagValue;
	this->makeInactive = false;
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

void SphericalBoundary::process(Candidate *candidate) const {
	double d = (candidate->current.getPosition() - center).mag();
	if (d >= radius) {
		candidate->setProperty(flag, flagValue);
		if (makeInactive) {
			candidate->setActive(false);
			candidate->setProperty("Deactivated", getDescription());
		}
	}
	if (limitStep)
		candidate->limitNextStep(radius - d + margin);
}

void SphericalBoundary::updateDescription() {
	std::stringstream s;
	s << "Spherical Boundary: radius " << radius << " around " << center;
	s << " Flag: " << flag << " -> " << flagValue;
	setDescription(s.str());
}

EllipsoidalBoundary::EllipsoidalBoundary(Vector3 focalPoint1,
		Vector3 focalPoint2, double majorAxis, std::string flag,
		std::string flagValue) {
	this->focalPoint1 = focalPoint1;
	this->focalPoint2 = focalPoint2;
	this->majorAxis = majorAxis;
	this->flag = flag;
	this->flagValue = flagValue;
	this->makeInactive = false;
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

void EllipsoidalBoundary::process(Candidate *candidate) const {
	Vector3 pos = candidate->current.getPosition();
	double d = (pos - focalPoint1).mag() + (pos - focalPoint2).mag();
	if (d >= majorAxis) {
		candidate->setProperty(flag, flagValue);
		if (makeInactive) {
			candidate->setActive(false);
			candidate->setProperty("Deactivated", getDescription());
		}
	}
	if (limitStep)
		candidate->limitNextStep(majorAxis - d + margin);
}

void EllipsoidalBoundary::updateDescription() {
	std::stringstream s;
	s << "Ellipsoidal Boundary: F1 = " << focalPoint1 / Mpc << ", F2 = "
			<< focalPoint2 / Mpc << ", major axis = " << majorAxis / Mpc;
	s << " Flag: " << flag << " -> " << flagValue;
	setDescription(s.str());
}

} // namespace mpc

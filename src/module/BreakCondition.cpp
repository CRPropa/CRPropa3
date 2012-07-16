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

void MaximumTrajectoryLength::process(Candidate *candidate) const {
	double l = candidate->getTrajectoryLength();
	if (l >= maxLength) {
		candidate->setActive(false);
		candidate->setProperty("Deactivated", getDescription());
	} else {
		candidate->limitNextStep(maxLength - l);
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

void MinimumEnergy::process(Candidate *candidate) const {
	if (candidate->current.getEnergy() <= minEnergy) {
		candidate->setActive(false);
		candidate->setProperty("Deactivated", getDescription());
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

void SmallObserverSphere::process(Candidate *candidate) const {
	double d = (candidate->current.getPosition() - center).getMag();
	if (d <= radius * 1.01) {
		candidate->setProperty(flag, flagValue);
		if (makeInactive) {
			candidate->setActive(false);
			candidate->setProperty("Deactivated", getDescription());
		}
	}
	candidate->limitNextStep((d - radius));
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

void PeriodicBox::process(Candidate *candidate) const {
	Vector3d position = candidate->current.getPosition();
	Vector3d n = ((position - origin) / size).floor(); // integers for translation
	if ((n.x != 0) or (n.y != 0) or (n.z != 0)) {
		candidate->current.setPosition(position - n * size);
		candidate->initial.setPosition(
				candidate->initial.getPosition() - n * size);
	}
}

void PeriodicBox::updateDescription() {
	std::stringstream s;
	s << "Periodic box: origin " << origin << ", size " << size;
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

void CubicBoundary::process(Candidate *candidate) const {
	Vector3d r = candidate->current.getPosition() - origin;
	double lo = r.min();
	double hi = r.max();
	if ((lo <= 0) or (hi >= size)) {
		candidate->setProperty(flag, flagValue);
		if (makeInactive) {
			candidate->setActive(false);
			candidate->setProperty("Deactivated", getDescription());
		}
	} else if (limitStep) {
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

void SphericalBoundary::process(Candidate *candidate) const {
	double d = (candidate->current.getPosition() - center).getMag();
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

void EllipsoidalBoundary::process(Candidate *candidate) const {
	Vector3d pos = candidate->current.getPosition();
	double d = pos.getDistanceTo(focalPoint1) + pos.getDistanceTo(focalPoint2);
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

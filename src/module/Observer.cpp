#include "mpc/module/Observer.h"

#include <sstream>

namespace mpc {

SmallObserverSphere::SmallObserverSphere(Vector3d center, double radius,
		std::string flag, std::string flagValue, bool makeInactive) :
		center(center), radius(radius), flag(flag), flagValue(flagValue), makeInactive(
				makeInactive) {
}

void SmallObserverSphere::process(Candidate *candidate) const {
	// current distance to observer sphere center
	double d = (candidate->current.getPosition() - center).getMag();

	// conservatively limit next step to prevent overshooting
	candidate->limitNextStep(fabs(d - radius));

	// no detection if outside of observer sphere
	if (d > radius)
		return;

	// previous distance to observer sphere center
	double dprev = (candidate->previous.getPosition() - center).getMag();

	// if particle was inside of sphere in previous step it has already been detected
	if (dprev <= radius)
		return;

	// else: detection
	candidate->setProperty(flag, flagValue);
	if (makeInactive)
		candidate->setActive(false);
}

void SmallObserverSphere::setCenter(Vector3d c) {
	center = c;
}

void SmallObserverSphere::setRadius(double r) {
	radius = r;
}

void SmallObserverSphere::setFlag(std::string f, std::string v) {
	flag = f;
	flagValue = v;
}

void SmallObserverSphere::setMakeInactive(bool b) {
	makeInactive = b;
}

std::string SmallObserverSphere::getDescription() const {
	std::stringstream s;
	s << "Small observer sphere: " << radius / Mpc;
	s << " Mpc radius around " << center / Mpc;
	s << " Mpc, Flag: '" << flag << "' -> '" << flagValue << "'";
	if (makeInactive)
		s << ", render inactivate";
	return s.str();
}

LargeObserverSphere::LargeObserverSphere(Vector3d center, double radius,
		std::string flag, std::string flagValue, bool makeInactive) :
		center(center), radius(radius), flag(flag), flagValue(flagValue), makeInactive(
				makeInactive) {
}

void LargeObserverSphere::process(Candidate *candidate) const {
	// current distance to observer sphere center
	double d = (candidate->current.getPosition() - center).getMag();

	// conservatively limit next step size to prevent overshooting
	candidate->limitNextStep(fabs(radius - d));

	// no detection if inside observer sphere
	if (d < radius)
		return;

	// previous distance to observer sphere center
	double dprev = (candidate->previous.getPosition() - center).getMag();

	// if particle was outside of sphere in previous step it has already been detected
	if (dprev >= radius)
		return;

	// else: detection
	candidate->setProperty(flag, flagValue);
	if (makeInactive)
		candidate->setActive(false);
}

void LargeObserverSphere::setCenter(Vector3d c) {
	center = c;
}

void LargeObserverSphere::setRadius(double r) {
	radius = r;
}

void LargeObserverSphere::setFlag(std::string f, std::string v) {
	flag = f;
	flagValue = v;
}

void LargeObserverSphere::setMakeInactive(bool b) {
	makeInactive = b;
}

std::string LargeObserverSphere::getDescription() const {
	std::stringstream s;
	s << "Large observer sphere: " << radius / Mpc;
	s << " Mpc radius around " << center / Mpc;
	s << " Mpc, Flag: '" << flag << "' -> '" << flagValue << "'";
	if (makeInactive)
		s << ", render inactivates";
	return s.str();
}

Observer1D::Observer1D() {
	setDescription("1D observer");
}

void Observer1D::process(Candidate *candidate) const {
	double x = candidate->current.getPosition().x;
	if (x > 0) {
		candidate->limitNextStep(x);
		return;
	}
	// else: detection
	candidate->setProperty("Detected", "");
	candidate->setActive(false);
}

DetectAll::DetectAll(std::string f, std::string v, bool m) :
		flag(f), flagValue(v), makeInactive(m) {
}

void DetectAll::process(Candidate *candidate) const {
	candidate->setProperty(flag, flagValue);
	if (makeInactive)
		candidate->setActive(false);
}

std::string DetectAll::getDescription() const {
	std::stringstream s;
	s << "DetectAll: Flag: " << flag << " -> " << flagValue;
	if (makeInactive)
		s << ", render inactive";
	return s.str();
}

} // namespace mpc

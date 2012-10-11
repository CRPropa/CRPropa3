#include "mpc/module/Observer.h"

#include <sstream>

namespace mpc {

SmallObserverSphere::SmallObserverSphere() :
		center(Vector3d(0, 0, 0)), radius(0), flag("Detected"), flagValue(""), makeInactive(
				true) {
}

SmallObserverSphere::SmallObserverSphere(Vector3d c, double r, std::string f,
		std::string v, bool b) :
		center(c), radius(r), flag(f), flagValue(v), makeInactive(b) {
}

void SmallObserverSphere::process(Candidate *c) const {
	double d = (c->current.getPosition() - center).getMag();
	if (d <= radius) {
		double dprev = (c->previous.getPosition() - center).getMag();
		if (dprev > radius) {
			c->setProperty(flag, flagValue);
			if (makeInactive)
				c->setActive(false);
		}
	}
	c->limitNextStep(fabs(d - radius));
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
		s << ", Inactivate";
	return s.str();
}

LargeObserverSphere::LargeObserverSphere() :
		center(Vector3d(0, 0, 0)), radius(0), flag("Detected"), flagValue(""), makeInactive(
				true) {
}

LargeObserverSphere::LargeObserverSphere(Vector3d c, double r, std::string f,
		std::string v, bool b) :
		center(c), radius(r), flag(f), flagValue(v), makeInactive(b) {
}

void LargeObserverSphere::process(Candidate *c) const {
	double d = (c->current.getPosition() - center).getMag();
	if (d >= radius) {
		double dprev = (c->previous.getPosition() - center).getMag();
		if (dprev < radius) {
			c->setProperty(flag, flagValue);
			if (makeInactive)
				c->setActive(false);
		}
	}
	c->limitNextStep(fabs(radius - d));
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
		s << ", Inactivate";
	return s.str();
}

Observer1D::Observer1D() {
	setDescription("1D observer");
}

void Observer1D::process(Candidate *c) const {
	double x = c->current.getPosition().x;
	if (x > 0) {
		c->limitNextStep(x);
		return;
	}
	// else: detection
	c->setProperty("Detected", "");
	c->setActive(false);
}

} // namespace mpc

#include "mpc/module/Observer.h"

#include <sstream>

namespace mpc {

SmallObserverSphere::SmallObserverSphere(Vector3d c, double r,
		std::string f, std::string v, bool b) :
		center(c), radius(r), flag(f), flagValue(v), makeInactive(b) {
	updateDescription();
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

void SmallObserverSphere::updateDescription() {
	std::stringstream s;
	s << "Small observer sphere: " << radius / Mpc;
	s << " Mpc radius around " << center / Mpc;
	s << " Mpc, Flag: '" << flag << "' -> '" << flagValue << "'";
	setDescription(s.str());
}

LargeObserverSphere::LargeObserverSphere(Vector3d c, double r,
		std::string f, std::string v, bool b) :
		center(c), radius(r), flag(f), flagValue(v), makeInactive(b) {
	updateDescription();
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

void LargeObserverSphere::updateDescription() {
	std::stringstream s;
	s << "Large observer sphere: " << radius / Mpc;
	s << " Mpc radius around " << center / Mpc;
	s << " Mpc, Flag: '" << flag << "' -> '" << flagValue << "'";
	setDescription(s.str());
}

} // namespace mpc

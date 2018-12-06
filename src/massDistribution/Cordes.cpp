#include "crpropa/massDistribution/Cordes.h"

#include <sstream>

namespace crpropa {

double Cordes::getHIIDensity(const Vector3d &position) const {
	double n=0;  // in m^-3

	double z=position.z;
	double R = sqrt(position.x*position.x+position.y*position.y);  //radius in galactic disk

	n += 0.025/ccm*exp(-fabs(z)/(1*kpc))*exp(-pow(R/(20*kpc),2));  // galactocentric component
	n += 0.2/ccm*exp(-fabs(z)/(0.15*kpc))*exp(-pow((R-4*kpc)/(2*kpc),2));  // anular component

	return n;
}

double Cordes::getDensity(const Vector3d &position) const {
	return Cordes::getHIIDensity(position);
}

double Cordes::getNucleonDensity(const Vector3d &position) const {
	return getHIIDensity(position);
}

bool Cordes::getIsForHI() {
	return isforHI;
}

bool Cordes::getIsForHII() {
	return isforHII;
}

bool Cordes::getIsForH2() {
	return isforH2;
}

std::string Cordes::getDescription() {
	std::stringstream s;
	s << "Density Cordes include HII component";
	return s.str();
}

}  // namespace crpropa

#include "crpropa/massDistribution/Massdistribution.h"

namespace crpropa {

void DensityList::addDensity(ref_ptr<Density> dens) {
	DensityList.push_back(dens);
}

double DensityList::getDensity(const Vector3d &position) const {
	double n = 0.;
	for (int i = 0; i < DensityList.size(); i++)
		n += DensityList[i]->getDensity(position);
	return n;
}

double DensityList::getHIDensity(const Vector3d &position) const {
	double n = 0.;
	for (int i = 0; i < DensityList.size(); i++)
		n += DensityList[i]->getHIDensity(position);
	return n;
}

double DensityList::getHIIDensity(const Vector3d &position) const {
	double n = 0.;
	for (int i = 0; i < DensityList.size(); i++)
		n += DensityList[i]->getHIIDensity(position);
	return n;
}

double DensityList::getH2Density(const Vector3d &position) const {
	double n = 0.;
	for (int i = 0; i < DensityList.size(); i++)
		n += DensityList[i]->getH2Density(position);
	return n;
}

double DensityList::getNucleonDensity(const Vector3d &position) const {
	double n = 0.;
	for (int i = 0; i < DensityList.size(); i++)
		n += DensityList[i]->getNucleonDensity(position);
	return n;
}

} //namespace crpropa


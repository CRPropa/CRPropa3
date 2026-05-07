#include "crpropa/massDistribution/Massdistribution.h"
#include <sstream>
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

std::string DensityList::getDescription() {
	std::stringstream ss; 
	ss << "DensityList with " << DensityList.size() << " modules: \n";
	for (int i = 0; i < DensityList.size(); i++) {
		ss << "density " << i + 1 << ": " << DensityList[i] -> getDescription();
	}
	
	return ss.str();
}

// ----------- DensityGrid -----------------------------------------------------------------

DensityGrid::DensityGrid(ref_ptr<Grid1f> grid, bool isForHI, bool isForHII, bool isForH2) : 
	grid(grid), isForHI(isForHI), isForHII(isForHII), isForH2(isForH2) {
		checkAndWarn();
	}

void DensityGrid::checkAndWarn() {
	bool allDeactivated = (isForHI == false) && (isForHII == false) && (isForH2 == false);
	if (allDeactivated) {
		KISS_LOG_WARNING << "DensityGrid has all types deactivated."
			<< "In this case all output will be n = 0. \n"
		 	<< "Please activate the intended particle type. \n";
	}
}

double DensityGrid::getHIDensity(const Vector3d &position) const {
	if (isForHI)
		return grid -> interpolate(position);
	else 
		return 0.;
}

double DensityGrid::getHIIDensity(const Vector3d &position) const {
	if (isForHII) 
		return grid -> interpolate(position);
	else
		return 0.;
}

double DensityGrid::getH2Density(const Vector3d &position) const {
	if (isForH2)
		return grid -> interpolate(position);
	else
		return 0.;
}

double DensityGrid::getDensity(const Vector3d &position) const {
	double n = 0;
	n += getHIDensity(position);
	n += getHIIDensity(position);
	n += getH2Density(position);

	return n;
}

double DensityGrid::getNucleonDensity(const Vector3d &position) const {
	double n = 0;
	n += getHIDensity(position);
	n += getHIIDensity(position);
	n += 2 * getH2Density(position);

	return n;
}

bool DensityGrid::getIsForHI() {
	return isForHI;
}

bool DensityGrid::getIsForHII() {
	return isForHII;
}

bool DensityGrid::getIsForH2() {
	return isForH2;
}

void DensityGrid::setIsForHI(bool b) {
	isForHI = b;
	checkAndWarn();
}

void DensityGrid::setIsForHII(bool b) {
	isForHII = b;
	checkAndWarn();
}

void DensityGrid::setIsForH2(bool b) {
	isForH2 = b;
	checkAndWarn();
}

void DensityGrid::setGrid(ref_ptr<Grid1f> grid) {
	this->grid = grid;
}

std::string DensityGrid::getDescription() {
	std::stringstream ss;
	ss << "Density in a given grid \n"; 
	return ss.str();
}

} //namespace crpropa
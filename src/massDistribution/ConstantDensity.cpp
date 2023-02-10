#include "crpropa/massDistribution/ConstantDensity.h"

#include "kiss/logger.h"

#include <sstream>

namespace crpropa{

ConstantDensity::ConstantDensity(double HI, double HII, double H2) {
	// set all types active which are not equal 0 and change number density
	if(HI!=0)
		setHI(true, HI);
	if(HII!=0)
		setHII(true, HII);
	if(H2!=0)
		setH2(true, H2);
}

double ConstantDensity::getDensity(const Vector3d &position) const {
	double n = 0;

	if(isHI)
		n += HIdensitynumber;
	if(isHII)
		n += HIIdensitynumber;
	if(isH2)
		n += H2densitynumber;

	// check if all densities are deactivated and raise warning if so
	if((isHI || isHII || isH2) == false){
		KISS_LOG_WARNING
			<< "\nCalled getNucleonDensity on fully deactivated ConstantDensity "
			<< "gas density model. In this case the density is allways set to 0. \n";
	}

	return n;
}

double ConstantDensity::getNucleonDensity(const Vector3d &position) const {
	double n = 0;

	if(isHI)
		n += HIdensitynumber;
	if(isHII)
		n += HIIdensitynumber;
	if(isH2)
		n += 2*H2densitynumber;

	// check if all densities are deactivated and raise warning if so
	if((isHI || isHII || isH2) == false){
		KISS_LOG_WARNING
			<< "\nCalled getNucleonDensity on fully deactivated ConstantDensity "
			<< "gas density model. In this case the density is allways set to 0. \n";
	}
	return n;
}

double ConstantDensity::getHIDensity(const Vector3d &position) const {
	return HIdensitynumber;
}

double ConstantDensity::getHIIDensity(const Vector3d &position) const{
	return HIIdensitynumber;
}

double ConstantDensity::getH2Density(const Vector3d &position) const{
	return H2densitynumber;
}

bool ConstantDensity::getIsForHI() {
	return isHI;
}

bool ConstantDensity::getIsForHII() {
	return isHII;
}

bool ConstantDensity::getIsForH2() {
	return isH2;
}

void ConstantDensity::setHI(bool activate, double densitynumber) {
	isHI = activate;
	HIdensitynumber = densitynumber;
}

void ConstantDensity::setHI(bool activate) {
	setHI(activate, HIdensitynumber);
}

void ConstantDensity::setHI(double densitynumber) {
	setHI(isHI, densitynumber);
}

void ConstantDensity::setHII(bool activate, double densitynumber) {
	isHII = activate;
	HIIdensitynumber = densitynumber;
}

void ConstantDensity::setHII(bool activate) {
	setHII(activate, HIIdensitynumber);
}

void ConstantDensity::setHII(double densitynumber) {
	setHII(isHII, densitynumber);
}

void ConstantDensity::setH2(bool activate, double densitynumber) {
	isH2 = activate;
	H2densitynumber = densitynumber;
}

void ConstantDensity::setH2(bool activate) {
	setH2(activate, H2densitynumber);
}

void ConstantDensity::setH2(double densitynumber) {
	setH2(isH2, densitynumber);
}

std::string ConstantDensity::getDescription() {
	std::stringstream s;
	s << "ConstantDensity:\n";
	s<< "HI component is ";
	if(!isHI)
		s<< "not ";
	s<< "active and has a density of " << HIdensitynumber/ccm << " cm^-3" << "\nHII component is ";
	if(!isHII)
		s<< "not ";
	s<<"active and has a density of " << HIIdensitynumber/ccm<<" cm^-3" <<  "\nH2 component is ";
	if(!isH2)
		s<<"not ";
	s<<"active and has a density of " << H2densitynumber/ccm << " cm^-3";
	return s.str();
}

}  // namespace crpropa

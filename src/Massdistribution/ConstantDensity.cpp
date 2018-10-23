#include "crpropa/Massdistribution/ConstantDensity.h"

namespace crpropa{

ConstantDensity::ConstantDensity(double HI, double HII, double H2) {

	//set all types activ which are not equal 0 and change densitynumber
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
	
	//check if any density is activ and give warning if not
	bool anyDensityActive = isHI||isHII||isH2;

	if(anyDensityActive == false){
		KISS_LOG_WARNING
			<< "\n called getDensity on deactivated ConstantDensity \n"
			<< "returned 0 density\n"
			<< "please use setHx(true, n)\n";
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
	
	//check if any density is activ and give warning if not
	bool anyDensityActive = isHI||isHII||isH2;

	if(anyDensityActive == false){
		KISS_LOG_WARNING
			<< "\n called getNucleonDensity on deactivated ConstantDensity \n"
			<< "returned 0 density\n"
			<< "please use setHx(true, n)\n";
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

bool ConstantDensity::getisforHI() {

	return isHI;
}

bool ConstantDensity::getisforHII() {

	return isHII;
}

bool ConstantDensity::getisforH2() {

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
	s << "ConstantDensity:";
	s<< "HI component is ";
	if(!isHI)
		s<< "not ";
	s<< "activ and has a density of " << HIdensitynumber/ccm << " cm^-3" << "      HII component is ";
	if(!isHII)
		s<< "not ";
	s<<"activ and  has a density of " << HIIdensitynumber/ccm<<" cm^-3" <<  "      H2 component is ";
	if(!isH2)
		s<<"not "; 
	s<<"activ and  has a density of " << H2densitynumber/ccm << " cm^-3";
	return s.str();
}

}//namespace 

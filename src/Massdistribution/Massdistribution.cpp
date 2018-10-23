#include "crpropa/Massdistribution/Massdistribution.h"


namespace crpropa {

void CustomDensity::add(ref_ptr<Density> dens) { 
	
	// check which density tpye is activated in loading density. Just use this part!		
	bool HI = dens->getisforHI();	
	bool HII= dens->getisforHII();
	bool H2 = dens->getisforH2();	
	
	bool nothingToLoad = !(HI || HII || H2);
	
	if(nothingToLoad == true){
		KISS_LOG_WARNING
		<<"\n tryed to add density to CustomDensity although no density type is activated. \n  nothing is loaded!\n";
		return ;
	}
	
	if(HI == true){
		HIDist = dens;
		isforHI = true;
		HIisload = true;
	}
	
	if(HII == true){
		HIIDist = dens;
		isforHII = true;
		HIIisload = true;
	}
	
	if(H2 == true){
		H2Dist = dens;
		isforH2 = true;
		H2isload = true;
	}
	
	return;
	
}


double CustomDensity::getDensity(const Vector3d &position) const{
	
	double n=0.;

	// warning if nothing is load in 	
	bool nothingLoadIn = !(HIisload || HIIisload || H2isload);
	if(nothingLoadIn)
	{	
		KISS_LOG_WARNING
		<< "\n called getDensity in CustomDensity alltough no option is loaded in. \n"
			<< "returned 0 density\n"
			<< "please load density\n";
		return 0;
	} 
	
	
	if(isforHI){
		n = HIDist->getHIDensity(position);
	}
	if(isforHII){
		n += this->HIIDist->getHIIDensity(position);
	}
	if(isforH2){
		n += this->H2Dist->getH2Density(position);
	}
	
	// warning if no option is activ in
	bool active = isforHI || isforHII || isforH2;
	if(active == false){
		KISS_LOG_WARNING
			<< "\n called getDensity on deactivated CustomDensity \n"
			<< "returned 0 density\n"
			<< "please activate\n";
	}
	
	return n;
}

double CustomDensity::getNucleonDensity(const Vector3d &position) const{
	
	double n=0.;
	bool nothingLoadIn = !(HIisload || HIIisload || H2isload);
	
	if(nothingLoadIn)
	{	
		KISS_LOG_WARNING
		<< "\n called getNucleonDensity in CustomDensity alltough no option is loaded in. \n"
			<< "returned 0 density\n"
			<< "please load density\n";
	} 

	
	if(isforHI){
		n = HIDist->getHIDensity(position);
	}
	if(isforHII){
		n += this->HIIDist->getHIIDensity(position);
	}
	if(isforH2){
		n += 2*this->H2Dist->getH2Density(position);
	}
	
	// warning if no option is activ
	bool activ = isforHI || isforHII || isforH2;
	if(activ == false){
		KISS_LOG_WARNING
			<< "\n called getNucleonDensity on deactivated CustomDensity. \n"
			<< "returned 0 density\n"
			<< "please activated\n";
	}
	
	return n;
}

double CustomDensity::getHIDensity(const Vector3d &position) const {
	
	if(HIisload)
		return HIDist->getHIDensity(position);
	
	// warning if no HI is load
	KISS_LOG_WARNING 
	<< "\n tryed to get HIDensity in CustomDensity, allthough no HI option is load in. \n"
	<< "Please load HI Option or use constant Density of 0 \n"
	<< "return a density of 0.\n";
	return 0.;
}

double CustomDensity::getHIIDensity(const Vector3d &position) const {
	
	if(HIIisload)
		return HIIDist->getHIIDensity(position);
	
	// warning if no HII is load
	KISS_LOG_WARNING 
	<< "\n tryed to get HIIDensity in CustomDensity, allthough no HII option is load in. \n"
	<< "Please load HII Option or use constant Density of 0 \n"
	<< "return a density of 0.\n";
	return 0.;
}

double CustomDensity::getH2Density(const Vector3d &position) const {
		
	if(H2isload)
		return H2Dist->getH2Density(position);
	
	// warning if no H2 is load
	KISS_LOG_WARNING 
	<< "\n tryed to get H2Density in CustomDensity, allthough no H2 option is load in. \n"
	<< "Please load H2 Option or use constant Density of 0 \n"
	<< "return a density of 0.\n";
	return 0.;
}

bool CustomDensity::getisforHI() {
	return isforHI;
}

bool CustomDensity::getisforHII() {
	return isforHII;
}

bool CustomDensity::getisforH2() {
	return isforH2;
}

void CustomDensity::setisforHI(bool HI) {
	
	if(HIisload==false && HI == true)
	{	
		KISS_LOG_WARNING
		<<"\n CustomDensity tryed to set HI true although no model is loaded. \n";
		return;
	}
	isforHI=HI;
	return;
}

void CustomDensity::setisforHII(bool HII) {

	if(HIIisload==false && HII == true)
	{	
		KISS_LOG_WARNING
		<<"\n CustomDensity tryed to set HII true although no model is loaded. \n";
		return;
	}
	isforHII=HII;
	return;
}

void CustomDensity::setisforH2(bool H2) {
	
	if(H2isload==false && H2 == true)
	{	
		KISS_LOG_WARNING
		<<"\n CustomDensity tryed to set H2 true although no model is loaded. \n";
		return;
	}
	isforH2=H2;
	return;
}

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


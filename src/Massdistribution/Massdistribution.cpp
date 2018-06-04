#include "crpropa/Massdistribution/Massdistribution.h"


namespace crpropa {
	

double Massdistribution::getDensity(const Vector3d &position) const{
	
	double n=0.;
	if(isforHI){
		n += distributionList[0]->getHIDensity(position);
	}
	if(isforHII){
		n += distributionList[1]->getHIIDensity(position);
	}
	if(isforH2){
		n += distributionList[2]->getH2Density(position);
	}
	return n;
}

void Massdistribution::add(Density &density) { //erase actual distribution and adds new 	
	
	bool HI = density.getisforHI();
	bool HII = density.getisforHII();
	bool H2 = density.getisforH2();
	auto it = distributionList.begin();
	if(HI){
		
		distributionList.erase(it);	
		distributionList.insert(it, &density);
		isforHI = HI;
	}
	if(density.getisforHII()){
		distributionList.erase(it+1);
		distributionList.insert(it+1, &density);
		isforHII = HII;
	}
	if(density.getisforH2()){
		distributionList.erase(it+2);
		distributionList.insert(it+2, &density);
		isforH2 = H2;
	}
	return;
}


	

} //namespace crpropa


#include "crpropa/Massdistribution/Cordes.h"

namespace crpropa {

double Cordes::getDensity(const Vector3d &position) const {

	return getHIIDensity(position);
}

double Cordes::getHIIDensity(const Vector3d &position) const {

	double n=0;
	
	double x=position.x/kpc;
	double y=position.y/kpc;
	double z=position.z/kpc;
	
	double R = sqrt(pow(x,2)+pow(y,2));
	
	n += 0.025*exp(-fabs(z)/1)*exp(-pow(R/20,2));
	n += 0.2*exp(-fabs(z)/0.15)*exp(-pow(R/2,2));
	
	return n;
}






}

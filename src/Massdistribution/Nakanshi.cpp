#include "crpropa/Massdistribution/Nakanshi.h"

namespace crpropa{ 

Nakanshi::Nakanshi(){
	isforHI = true;
	isforHII = false;
	isforH2 = true;		//model typ: do not change!
}

double Nakanshi::getDensity(const Vector3d &position) {
	double n = 0;
	n += getHIDensity(position);
	n += getH2Density(position);
	return n;
}

double Nakanshi::getHIDensity(const Vector3d &position) {
	
	double z = position.z/kpc;
	double planedensity = getHIPlanedensity(position);
	double scaleheight = getHIScaleheight(position);
	return planedensity*pow(0.5,pow(z/scaleheight,2));
}

double Nakanshi::getHIScaleheight(const Vector3d &position) {
		
	double R = sqrt(pow(position.x,2)+pow(position.y,2));	//Galactic Radius		
	double scaleheight = 1.06*(116.3 +19.3*R+4.1*pow(R,2)-0.05*pow(R,3));
	return scaleheight;
	}

double Nakanshi::getHIPlanedensity(const Vector3d &position) {

	double R = sqrt(pow(position.x,2)+pow(position.y,2));
	double planedensity = 0.94*(0.6*exp(-R/2.4)+0.24*exp(-pow((R-9.5)/4.8,2)));
	return planedensity;
	}

double Nakanshi::getH2Density(const Vector3d &position) {
	
	double z = position.z/kpc;
	double plane = getH2Planedensity(position);
	double scaleheight = getH2Scaleheight(position);
	return 2*plane*pow(0.5,pow(position.z,2));	//Factor 2 for both nuclei
}
	
double Nakanshi::getH2Scaleheight(const Vector3d &position)  {

	double R = sqrt(pow(position.x,2)+ pow(position.y,2));
	double Scaleheight = 1.06e-3*( 11.2*exp(0.28*R)+42.78);
	return Scaleheight;
}

double Nakanshi::getH2Planedensity(const Vector3d &position) {

	double R = sqrt(pow(position.x,2)+pow(position.y,2));
	double Planedensity = 11.2*exp(-pow(R,2)/0.85) +0.83*exp(-pow((R-4)/3.2,2));
	return 0.94*Planedensity;
}

std::string Nakanshi::getDiscription() {
	std::string Description="modell Nakanshi ";
	if(isforHI)
	{
		Description += "2003 for HI";
		if(isforH2)
		{
			Description += "and Nakanshi 2006 for H2";
		}
	return Description;
	}
	if(isforH2)
	{
		Description += "2006 for H2";
	}
	else
	{
		Description += "not in use";
	}
	return Description;
}

} //namespace crpropa










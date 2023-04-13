#include "crpropa/massDistribution/Ferriere.h"
#include "crpropa/Common.h"

#include "kiss/logger.h"

#include <sstream>

namespace crpropa {

Vector3d Ferriere::CMZTransformation(const Vector3d &position) const {
	// set galactocentric coordinate system with the Sun at (-8.5,0.,0.) instead of (8.5, 0, 0) to be consistand with JF12 implementation
	double x = -position.x;
	double y = -position.y;

	double xC = -50*pc;		//offset
	double yC = 50*pc;
	double sinTc = sin(70.*deg);
	double cosTc = cos(70.*deg);

	Vector3d pos;
	pos.x = (x - xC)*cosTc + (y - yC)*sinTc;
	pos.y = -(x - xC)*sinTc + (y - yC)*cosTc;
	pos.z = position.z;

	return pos;
}

Vector3d Ferriere::DiskTransformation(const Vector3d &position) const {
	// set galactocentric coordinate system with the Sun at (-8.5,0.,0.) instead of (8.5, 0, 0) to be consistand with JF12 implementation
	double x = -position.x;
	double y = - position.y;
	double z = position.z;

	double alphaD = 13.5*deg;  // rotation arround x-axis
	double sinAd = sin(alphaD);
	double cosAd = cos(alphaD);
	double betaD = 20.*deg;  // rotation arround y'-axis
	double sinBd = sin(betaD);
	double cosBd = cos(betaD);
	double TettaD = 48.5*deg;  // rotation arround x"-axis
	double sinTd = sin(TettaD);
	double cosTd = cos(TettaD);

	Vector3d pos;

	pos.x = x*cosBd*cosTd - y*(sinAd*sinBd*cosTd -cosAd*sinTd)-z*(cosAd*sinBd*cosTd +sinAd*sinTd);

	pos.y =  -x*cosBd*sinTd;
	pos.y += y*(sinAd*sinBd*sinTd +cosAd*cosTd);
	pos.y += z*(cosAd*sinBd*sinTd -sinAd*cosTd);

	pos.z = x*sinBd;
	pos.z += y*sinAd*cosBd;
	pos.z += z*cosAd*cosBd;

	return pos;
}

double Ferriere::getHIDensity(const Vector3d &position) const {
	double n = 0;
	double R = sqrt(position.x*position.x+position.y*position.y);

	if(R<3*kpc)
	{
		// density at center
		Vector3d pos = CMZTransformation(position);  // coordinate trafo
		double x = pos.x/pc;  // all units in pc
		double y = pos.y/pc;
		double z = pos.z/pc;

		double A = sqrt(x*x+2.5*2.5*y*y);
		double nCMZ = 8.8/ccm*exp(-pow_integer<4>((A-125.)/137))*exp(-pow_integer<2>(z/54.));

		// density in disk
		pos = DiskTransformation(position);  // Corrdinate Trafo
		x = pos.x/pc;  // all units in pc
		y = pos.y/pc;
		z = pos.z/pc;

		A = sqrt(x*x+3.1*3.1*y*y);
		double nDisk = 0.34/ccm*exp(-pow_integer<4>((A-1200.)/438.))*exp(-pow_integer<2>(z/120));

		n = nCMZ + nDisk;
	}
	else{  // outer region
		double z = position.z/pc;
		double a;
		if(R<=Rsun){
			a = 1;
		} else {
			a = R/Rsun;
		}

		double nCold = 0.859*exp(-pow_integer<2>(z/(127*a))); // cold HI
		nCold += 0.047*exp(-pow_integer<2>(z/(318*a)));
		nCold += 0.094*exp(-fabs(z)/(403*a));
		nCold *= 0.340/ccm/(a*a);

		double nWarm = (1.745 - 1.289/a)*exp(-pow_integer<2>(z/(127*a)));  // warm HI
		nWarm += (0.473 - 0.070/a)*exp(-pow_integer<2>(z/(318*a)));
		nWarm += (0.283 - 0.142/a)*exp(-fabs(z)/(403*a));
		nWarm *= 0.226/ccm/a;

		n = nWarm + nCold;
	}

	return n;
}

double Ferriere::getHIIDensity(const Vector3d &position) const {
	double n = 0;
	double R = sqrt(position.x*position.x+position.y*position.y);

	if(R< 3*kpc){   // inner
		double x = position.x/pc;
		double y = position.y/pc;
		double z = position.z/pc;

		// warm interstellar matter
		double H = (x*x + pow_integer<2>(y+10.))/(145*145);
		double nWIM = exp(-H)* exp(-pow_integer<2>(z+20.)/(26*26));
		nWIM += 0.009*exp(-pow_integer<2>((R/pc-3700)/(0.5*3700)))*1/pow_integer<2>(cosh(z/140.));
		nWIM += 0.005*cos(M_PI*R/pc*0.5/17000)*1/pow_integer<2>(cosh(z/950.));
		nWIM *= 8.0/ccm;  // rescaling with density at center

		//very hot interstellar matter
		double alphaVH = 21.*deg;  // angle for very hot IM in radiant
		double cosA = cos(alphaVH);
		double sinA = sin(alphaVH);
		double etta = y*cosA+z*sinA;  // coordinate transformation for VHIM along major axis
		double chi = -y*sinA+z*cosA;

		double nVHIM = 0.29/ccm*exp(-((x*x+etta*etta)/(162.*162.)+chi*chi/(90*90)));

		n = nWIM + nVHIM;
	} else {  // outer region
		double z = position.z;

		double nWarm = 0.0237/ccm*exp(-(R*R-Rsun*Rsun)/pow_integer<2>(37*kpc))*exp(-fabs(z)/(1*kpc));
		nWarm += 0.0013/ccm*exp(-(pow_integer<2>(R-4*kpc)-pow_integer<2>(Rsun-4*kpc))/pow_integer<2>(2*kpc))*exp(-fabs(z)/(150*pc));

		double nHot = 0.12*exp(-(R-Rsun)/(4.9*kpc));
		nHot += 0.88*exp(-(pow_integer<2>(R-4.5*kpc)-pow_integer<2>(Rsun-4.5*kpc))/pow_integer<2>(2.9*kpc));
		nHot *= pow(R/Rsun, -1.65);
		nHot *= exp(-fabs(z)/((1500*pc)*pow(R/Rsun,1.65)));
		nHot *= 4.8e-4/ccm;

		n= nWarm + nHot;
	}
	return n;
}

double Ferriere::getH2Density(const Vector3d &position) const{
	double n=0;
	double R=sqrt(position.x*position.x+position.y*position.y);

	if(R<3*kpc) {
		// density at center
		Vector3d pos =CMZTransformation(position);  // coord trafo for CMZ
		double x = pos.x/pc;  // all units in pc
		double y = pos.y/pc;
		double z = pos.z/pc;

		double A = sqrt(x*x+pow(2.5*y,2));  // ellipticity
		double nCMZ = exp(-pow((A-125.)/137.,4))*exp(-pow(z/18.,2));
		nCMZ *= 150/ccm;  // rescaling

		// density in disk
		pos = DiskTransformation(position);  // coord trafo for disk
		x=pos.x/pc;
		y=pos.y/pc;
		z=pos.z/pc;

		A = sqrt(x*x+pow_integer<2>(3.1*y));
		double nDISK = exp(-pow_integer<4>((A-1200)/438))*exp(-pow_integer<2>(z/42));
		nDISK *= 4.8/ccm;  // rescaling

		n = nCMZ + nDISK;
	} else {  // outer region
		double z = position.z/pc;
		n = pow(R/Rsun, -0.58);
		n *= exp(-(pow_integer<2>(R-4.5*kpc)-pow_integer<2>(Rsun-4.5*kpc))/pow_integer<2>(2.9*kpc));
		n *= exp(-pow_integer<2>(z/(81*pow(R/Rsun,0.58))));
		n *= 0.58/ccm;  // rescaling
	}

	return n;
}

double Ferriere::getDensity(const Vector3d &position) const{
	double n=0;
	if(isforHI){
		n += getHIDensity(position);
	}
	if(isforHII){
		n+=getHIIDensity(position);
	}
	if(isforH2){
		n+=getH2Density(position);
	}
	// check if all densities are deactivated and raise warning if so
	if((isforHI || isforHII || isforH2) == false){
		KISS_LOG_WARNING
			<< "\nCalled getDensity on fully deactivated Ferriere \n"
			<< "gas density model. The total density is set to 0.";
	}

	return n;
}

double Ferriere::getNucleonDensity(const Vector3d &position) const{
	double n=0;
	if(isforHI){
		n += getHIDensity(position);
	}
	if(isforHII){
		n+=getHIIDensity(position);
	}
	if(isforH2){
		n+= 2*getH2Density(position);
	}

	// check if all densities are deactivated and raise warning if so
	if((isforHI || isforHII || isforH2) == false){
		KISS_LOG_WARNING
			<< "\nCalled getNucleonDensity on fully deactivated Ferriere \n"
			<< "gas density model. The total density is set to 0.";
	}

	return n;
}

void Ferriere::setIsForHI(bool HI){
	isforHI = HI;
}

void Ferriere::setIsForHII(bool HII){
	isforHII = HII;
}

void Ferriere::setIsForH2(bool H2){
	isforH2 = H2;
}

bool Ferriere::getIsForHI(){
	return isforHI;
}

bool Ferriere::getIsForHII(){
	return isforHII;
}

bool Ferriere::getIsForH2(){
	return isforH2;
}

std::string Ferriere::getDescription() {
	std::stringstream s;
	s << "Density model Ferriere 2007:\n";
	s<< "HI component is ";
	if(!isforHI)
		s<< "not ";
	s<< "active.\nHII component is ";
	if(!isforHII)
		s<< "not ";
	s<<"active.\nH2 component is ";
	if(!isforH2)
		s<<"not ";
	s<<"active.";
	return s.str();
}

}  // namespace crpropa
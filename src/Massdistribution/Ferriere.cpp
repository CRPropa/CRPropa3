#include "crpropa/Massdistribution/Ferriere.h"

namespace crpropa {


Vector3d Ferriere::CMZTrafo(const Vector3d &position) const {
	
	// set galactocentric coordinate system with the Sun at (-8.5,0.,0.) instead of (8.5, 0, 0) to be consistand with JF12 implementation
	double x = -position.x;
	double y = -position.y;
	
	double xC = -50*pc;		//offset
	double yC = 50*pc;	
	double sinTc = sin(70./180.*M_PI);	//Tc = 70 deg
	double cosTc = cos(70./180.*M_PI);
	
	
	Vector3d pos;
	pos.x = (x - xC)*cosTc + (y - yC)*sinTc;
	pos.y = -(x - xC)*sinTc + (y - yC)*cosTc;
	pos.z = position.z;
		
	
	return pos;
}

Vector3d Ferriere::DISKTrafo(const Vector3d &position) const { 

	// set galactocentric coordinate system with the Sun at (-8.5,0.,0.) instead of (8.5, 0, 0) to be consistand with JF12 implementation
	double x = -position.x;
	double y = - position.y;
	double z = position.z;
	
	double alphaD = 13.5/180.*M_PI;	// rotation arround x-axis
	double sinAd = sin(alphaD);
	double cosAd = cos(alphaD);
	double betaD = 20./180.*M_PI; // rotation arround y'-axis
	double sinBd = sin(betaD);
	double cosBd = cos(betaD);
	double TettaD = 48.5/180.*M_PI; // rotation arround x"-axis
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
	bool inner = R<3*kpc;
	
	if(inner == true) 
	{
		double nCMZ = 0;	//center
		
		Vector3d pos = CMZTrafo(position);	//coordinate trafo
		double x = pos.x/pc;	// all units in pc
		double y = pos.y/pc;
		double z = pos.z/pc;
		
		double A = sqrt(x*x+2.5*2.5*y*y);
		nCMZ = 8.8*exp(-pow((A-125.)/137,4))*exp(-pow(z/54.,2));
		
		double nDisk = 0;		//disk
		
		pos = DISKTrafo(position);	//Koordinaten Trafo
		x = pos.x/pc;	// all units in pc
		y = pos.y/pc;
		z = pos.z/pc;
		
		A = sqrt(x*x+3.1*3.1*y*y);
		nDisk = 0.34*exp(-pow((A-1200.)/438.,4))*exp(-pow(z/120,2));
		
		n = nCMZ + nDisk;
	}
	else{ // outer region

		double z = position.z/pc;	
		double a;
		if(R<=Rsun){
			a = 1;
		}
		else {
			a = R/Rsun;
		}
		
		
		double nCold = 0;	//cold HI
		nCold += 0.859*exp(-pow(z/(127*a),2));
		nCold += 0.047*exp(-pow(z/(318*a),2));
		nCold += 0.094*exp(-fabs(z)/(403*a));
		nCold *= 0.340/(a*a);
		
		
		double nWarm = 0;	//warm HI
		nWarm += (1.745 - 1.289/a)*exp(-pow(z/(127*a),2));
		nWarm += (0.473 - 0.070/a)*exp(-pow(z/(318*a),2));
		nWarm += (0.283 - 0.142/a)*exp(-fabs(z)/(403*a));
		nWarm *= 0.226/a;
		
		n = nWarm + nCold;

	}
	
	return n/ccm;
}

double Ferriere::getHIIDensity(const Vector3d &position) const {

	double n = 0;
	double R = sqrt(position.x*position.x+position.y*position.y);
	bool inner = R< 3*kpc;
		
	if(inner == true){   //inner
	
	double x = position.x/pc;
	double y = position.y/pc;
	double z = position.z/pc;
	
	//warm interstellar matter
	double nWIM = 0;		
	
	double H = (x*x + pow(y+10.,2))/(145*145);
	nWIM += exp(-H);
	double n0 = nWIM;
	nWIM *= exp(-pow(z+20.,2)/(26*26));
	double n1 = nWIM;
	nWIM += 0.009*exp(-pow((R/pc-3700)/(0.5*3700),2))*1/pow(cosh(z/140.),2);
	double n2 = nWIM;
	nWIM += 0.005*cos(M_PI*R/pc*0.5/17000)*1/pow(cosh(z/950.),2);
	double n3 = nWIM;
	
	nWIM *= 8.0;	//rescaling
	
	
	//very hot interstellar matter
	double nVHIM = 0;
	double alphaVH = 21./180*M_PI;		//angle for very hot IM in radiant 
	double cosA = cos(alphaVH);
	double sinA = sin(alphaVH);
	double etta = y*cosA+z*sinA;		// coordinate transformation for VHIM along major axis
	double chi = -y*sinA+z*cosA;
	
	nVHIM = 0.29*exp(-((x*x+etta*etta)/(162.*162.)+chi*chi/(90*90)));
	
	n = nWIM + nVHIM;

	}
	else {		// outer region
		
		double z = position.z;
		
		double nWarm = 0;
		nWarm += 0.0237*exp(-(pow(R,2)-pow(Rsun,2))/pow(37*kpc,2))*exp(-fabs(z)/(1*kpc));
		nWarm += 0.0013*exp(-(pow(R-4*kpc,2)-pow(Rsun-4*kpc,2))/pow(2*kpc,2))*exp(-fabs(z)/(150*pc));
		
		double nHot = 0;
		nHot += 0.12*exp(-(R-Rsun)/(4.9*kpc));
		nHot += 0.88*exp(-(pow(R-4.5*kpc,2)-pow(Rsun-4.5*kpc,2))/pow(2.9*kpc,2));
		nHot *= pow(R/Rsun, -1.65);
		nHot *= exp(-fabs(z)/((1500*pc)*pow(R/Rsun,1.65)));
		nHot *= 4.8e-4;
		
		n= nWarm + nHot;
		
	}
		
	return n/ccm;
}


double Ferriere::getH2Density(const Vector3d &position) const{

	double n=0;
	double R=sqrt(pow(position.x,2)+pow(position.y,2));
	bool inner = R<3*kpc;	
	
	if(inner == true) {		
	
		double nCMZ = 0;
		double nDISK = 0;
		
		Vector3d pos =CMZTrafo(position); // coord trafo for CMZ
		double x = pos.x/pc;	// all units in pc
		double y = pos.y/pc;
		double z = pos.z/pc;
		
		double A = sqrt(pow(x,2)+pow(2.5*y,2));	// ellipticity
		nCMZ = exp(-pow((A-125.)/137.,4))*exp(-pow(z/18.,2));
		nCMZ *= 150;		// rescaling
		
		pos = DISKTrafo(position);	// coord trafo for disk
		x=pos.x/pc;
		y=pos.y/pc;
		z=pos.z/pc;
		
		A = sqrt(pow(x,2)+pow(3.1*y,2));
		nDISK = exp(-pow((A-1200)/438,4))*exp(-pow(z/42,2));
		nDISK *= 4.8;		//rescaling
		
		n = nCMZ + nDISK;
	}	
	else {		//outer region
		double z = position.z/pc;
		n = pow(R/Rsun, -0.58);
		n *= exp(-(pow(R-4.5*kpc,2)-pow(Rsun-4.5*kpc,2))/pow(2.9*kpc,2));
		n *= exp(-pow(z/(81*pow(R/Rsun,0.58)),2));
		n *= 0.58;		//rescaling
	}
	
	return n/ccm;
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
	
	//check if any density is activ and give warning if not
	bool anyDensityActive = isforHI||isforHII||isforH2;

	if(anyDensityActive == false){
		KISS_LOG_WARNING
			<< "\n called getDensity on deactivated Ferriere \n"
			<< "returned 0 density\n"
			<< "please activate \n";
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
	
	//check if any density is activ and give warning if not
	bool anyDensityActive = isforHI||isforHII||isforH2;

	if(anyDensityActive == false){
		KISS_LOG_WARNING
			<< "\n called getDensity on deactivated ConstantDensity \n"
			<< "returned 0 density\n"
			<< "please activate\n";
	}
	
	return n;
}


void Ferriere::setisforHI(bool HI){

	isforHI = HI;
}

void Ferriere::setisforHII(bool HII){
	
	isforHII = HII;
}

void Ferriere::setisforH2(bool H2){
	
	isforH2 = H2;
}

bool Ferriere::getisforHI(){
	
	return isforHI;
}

bool Ferriere::getisforHII(){
	
	return isforHII;
}
	
bool Ferriere::getisforH2(){
	
	return isforH2;
}

std::string Ferriere::getDescription() {
	
	std::stringstream s;
	s << "Density modell Ferriere 2007: ";
	s<< "HI component is ";
	if(!isforHI)
		s<< "not ";
	s<< "activ. HII component is ";
	if(!isforHII)
		s<< "not ";
	s<<"activ. H2 component is ";
	if(!isforH2)
		s<<"not "; 
	s<<"activ.";
	return s.str();
}


} //namespace 

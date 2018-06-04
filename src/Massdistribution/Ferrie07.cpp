#include "crpropa/Massdistribution/Ferrie07.h"

namespace crpropa {

     Ferrie::Ferrie(){}

Vector3d Ferrie::CMZTrafo(const Vector3d &position) const {
	
	double xC = 50;		//offset
	double yC = 50;
	double ThettaC = 70/180*M_PI; //radiant
	
	Vector3d pos;
	pos.x = (position.x - xC)*cos(ThettaC) + (position.y -yC)*sin(ThettaC);
	pos.y = -(position.x -xC)*sin(ThettaC) + (position.y - yC)*cos(ThettaC);
	pos.z = position.z;
	return pos;
}

Vector3d Ferrie::DISKTrafo(const Vector3d &position) const { 

	double alphaD = 13.5/180*M_PI;
	double betaD = 20/180*M_PI;
	double TettaD = 48.5/180*M_PI;
	
	Vector3d pos;
	
	double x = position.x;
	double y = position.y;
	double z = position.z;
	
	pos.x = x*cos(betaD)*cos(TettaD) - y*(sin(alphaD)*sin(betaD)*cos(TettaD) -cos(alphaD)*sin(TettaD))-z*(cos(alphaD)*sin(betaD)*cos(TettaD) +sin(alphaD)*sin(TettaD));
	
	pos.y -= x*cos(betaD)*sin(TettaD);
	pos.y += y*(sin(alphaD)*sin(betaD)*sin(TettaD) +cos(alphaD)*cos(TettaD));
	pos.y += z*(cos(alphaD)*sin(betaD)*sin(TettaD) -sin(alphaD)*sin(TettaD));
	
	pos.z += x*sin(betaD);
	pos.z += y*sin(alphaD)*cos(betaD);
	pos.z += z*cos(alphaD)*cos(betaD);
	
	return pos;
}
	
	

    double Ferrie::getDensity( const crpropa::Vector3d &position) const{

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
	return n;
}

double Ferrie::getHIDensity(const Vector3d &position) const {
	
	double n = 0;
	double R = sqrt(pow(position.x/pc,2)+pow(position.y/pc,2));
	
	if(R<3*kpc) //INNEN
	{
		double nCMZ = 0;	//Center
		
		Vector3d pos = CMZTrafo(position);	//Koordinaten Trafo
		double x = pos.x;
		double y = pos.y;
		double z = pos.z;
		
		double A = sqrt(pow(x,2)+pow(2.5*y,2));
		nCMZ = 8.8*exp(-pow((A-250/2)/137,4))*exp(-pow(z/54,2));
		
		double nDisk = 0;		//Disk
		
		pos = DISKTrafo(position);	//Koordinaten Trafo
		x = pos.x/pc;
		y = pos.y/pc;
		z = pos.z/pc;
		
		A = sqrt(pow(x,2)+pow(3.1*y,2));
		nDisk = 0.34*exp(-pow((A-1200)/438,4))*exp(-pow(z/120,2));
		
		
		n = nCMZ + nDisk;
		
	}
	else{ //AUßEN

		double z = position.z/pc;	
		double a;
		if(R<=Rsun){
			a= 1;
		}
		else {
			a = R/Rsun;
		}
		
		
		double nCold =0;		//cold HI
		nCold += 0.859*exp(-pow(z/(127*pc*a),2));
		nCold += 0.047*exp(-pow(z/(318*pc*a),2));
		nCold += 0.094*exp(-fabs(z)/(403*pc*a));
		nCold *= 0.304/(pow(a,2));
		
		
		double nWarm =0;		//warm HI
		nWarm += (1.745 - 0.283/a)*exp(-pow(z/(127*pc*a),2));
		nWarm += (0.473 - 0.070/a)*exp(-pow(z/(318*pc*a),2));
		nWarm += (0.283 - 0.142/a)*exp(-fabs(z)/(403*pc*a));
		nWarm *= 0.226/a;
		
		
		n = nWarm + nCold;
	}
	return n;
}

double Ferrie::getHIIDensity(const Vector3d &position) const {

	double n = 0;
	double R = sqrt(pow(position.x/pc,2)+pow(position.y/pc, 2));
	if(R<3000*pc){   //innen
	
	double x = position.x/pc;
	double y = position.y/pc;
	double z = position.z/pc;
	
	double nWIM =0;		//warm IM
	double nVHIM = 0;	//very hot IM
	
	nWIM += exp(-(pow(x,2)+pow(y+10,2))/pow(145,2))*exp(-pow((z+20)/26,2));
	nWIM += 0.009*exp(-pow((R-3700)/(0.5*3750),2))*1/pow(cosh(z/140),2);
	nWIM += 0.005*cos(M_PI*R*0.5/17000)*1/pow(cosh(z/950),2);
	
	nWIM *= 8.0;	//center Density for scaling 
	
	double alphaVH = 21/180*M_PI;		//angel for very hot IM in radiant 
	double etta = y*cos(alphaVH)+z*sin(alphaVH);
	double chi = -y*sin(alphaVH)+z*cos(alphaVH);
	
	nVHIM = 0.29*exp(-(pow(x,2)+pow(etta,2))/pow(162,2)+pow(chi/90,2));
	
	n = nWIM + nVHIM;
	}
	else {		//außen
		
		double z = position.z/pc;
		
		double nWarm = 0;
		double nHot = 0;
		
		nWarm += 0.0237*exp(-(pow(R,2)-pow(Rsun,2))/pow(37000,2))*exp(-fabs(z)/1000);
		nWarm += 0.0013*exp(-(pow(R-4000,2)-pow(Rsun-4000,2))/pow(2000,2))*exp(-fabs(z)/150);
		
		nHot += 0.12*exp(-(R-Rsun)/4900);
		nHot += 0.88*exp(-(pow(R-4500,2)-pow(Rsun-4500,2))/pow(2900,2));
		
		nHot *= pow(R/Rsun, -1.65);
		nHot *= exp(-fabs(z)/(1500*pow(R/Rsun,1.65)));
		
		n= nWarm + nHot;
		
	}
	
	return n;
}


double Ferrie::getH2Density(const Vector3d &position) const{

	double n=0;
	double R=sqrt(pow(position.x,2)+pow(position.y,2));
	
	if(R<3000*pc) {		//innen
	
		double nCMZ = 0;
		double nDISK = 0;
		
		Vector3d pos = CMZTrafo(position);
		double x = pos.x/pc;
		double y = pos.y/pc;
		double z = pos.z/pc;
		
		double A = sqrt(pow(x,2)+pow(2.5*y,2));
		nCMZ = exp(-pow((A-125)/137,4))*exp(-fabs(z)/18);
		nCMZ *= 150;		//Density at Center for scale
		
		pos = DISKTrafo(position);
		x=pos.x/pc;
		y=pos.y/pc;
		z=pos.z/pc;
		
		A = sqrt(pow(x,2)+pow(3.1*y,2));
		nDISK = exp(-pow((A-1200)/438,4))*exp(-pow(z/42,2));
		nDISK *= 4.8;		//Density at Center for scale
		
		n = nCMZ + nDISK;
	}	
	else {		//außen
		
		double z = position.z/pc;
		n = pow(R/Rsun, -0.58);
		n *= exp(-(pow(R-4000,2)-pow(Rsun-4000,2))/pow(2900,2));
		n *= exp(-pow(z/(81*pow(R/Rsun,0.58)),2));
		n *= 0.58;		//Density at Center for scale
	}
	return n;
}

void Ferrie::setisforHI(bool HI){

	isforHI = HI;
	return;
}

void Ferrie::setisforHII(bool HII){
	
	isforHII = HII;
	return;
}

void Ferrie::setisforH2(bool H2){
	
	isforH2 = H2;
	return;
}

bool Ferrie::getisforHI(){
	
	return isforHI;
}

bool Ferrie::getisforHII(){
	
	return isforHII;
}
	
bool Ferrie::getisforH2(){
	
	return isforH2;
}

} //namespace 

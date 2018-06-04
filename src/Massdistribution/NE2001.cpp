#include "crpropa/Massdistribution/NE2001.h"

#include <fstream>
#include <algorithm> 


namespace crpropa {

NE2001::NE2001() {
	
	Rsun = 8.5*kpc;
	n1 = 0.033*pow(cm,-3);
	H1 = 0.95*kpc;
	A1 = 17*kpc;
	n2 = 0.09*pow(cm,-3);
	H2 = 0.14*kpc;
	A2 = 3.7*kpc;
	Aa = 8.5*kpc;	
		
	fjna[0] = 0.30*0.5*pow(cm,-3);
	fjna[1] = 0.3*1.2*pow(cm,-3);
	fjna[2] = 0.3*1.3*pow(cm,-3);
	fjna[3] = 0.3*1.0*pow(cm,-3);
	fjna[4] = 0.3*0.25*pow(cm,-3);
	
	hjHa[0] = 0.25*1.0*kpc;
	hjHa[1] = 0.25*0.8*kpc;
	hjHa[2] = 0.25*1.3*kpc;
	hjHa[3] = 0.25*1.5*kpc;
	hjHa[4] = 0.25*1.0*kpc;
	
	wjwa[0] = 0.6*1*kpc;
	wjwa[1] = 0.6*1.5*kpc;
	wjwa[2] = 0.6*1*kpc;
	wjwa[3] = 0.6*0.8*kpc;
	wjwa[4] = 0.6*1*kpc;
	
	nGC = 10.0*pow(cm,-3);
	RGC = 0.145*kpc;
	hGC = 0.026*kpc;
	PositionGC = Vector3d(-0.01*kpc, 0, -0.02*kpc);
	
	std::ifstream fin("share/crpropa/spiralArms.txt");
	
	if(!fin) {
		throw std::runtime_error("NE2001: SpiralArm data not found");
	}
	
	while(fin.peek() == '#')
		fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	
	for (int i=0;i<1000; i++) {
		for (int j = 0;j<5;j++) {
			
			fin >> xArm[i][j] >> yArm[i][j];
			if(fin.eof())
				throw std::runtime_error("load SpiralArm data: file too short");
		}
	}
	fin.close();
}

bool NE2001::getisforHI() {
	return isforHI;
}
bool NE2001::getisforHII() {
	return isforHII;
}
bool NE2001::getisforH2() {
	return isforH2;
}

bool NE2001::getuseSpiralArm() {
	return useSpiralArm;
}
bool NE2001::getuseThikDisk() {
	return useThikDisk;
}
bool NE2001::getuseThinDisk() {
	return useThinDisk;
}
bool NE2001::getuseGalacticCenter() {
	return useGalacticCenter;
}
/*
bool NE2001::getuseVoids() {
	return useVoids;
}
bool NE2001::getuseClumps() {
	return useClumps;
}
bool NE2001::getuseLISM() {
	return useLISM;
} */

void NE2001::setuseSpiralArm(bool use) {
	useSpiralArm = use;
}
void NE2001::setuseThikDisk(bool use) {
	useThikDisk = use;
}
void NE2001::setuseThinDisk(bool use) {
	useThinDisk = use;
}
void NE2001::setuseGalacticCenter(bool use) {
	useGalacticCenter = use;
}
/*
void NE2001::setuseVoids(bool use) {
	useVoids = use;
}
void NE2001::setuseClumps(bool use) {
	useClumps = use;
}
void NE2001::setuseLISM(bool use) {
	useLISM = use;
} */

double NE2001::getDensity(const Vector3d &position) const {
	return getHIIDensity(position);
}
	
double NE2001::getHIIDensity(const Vector3d &position) const {
	double n = 0;
	double x = position.x;
	double y = position.y;
	double r = sqrt(pow(x,2) +pow(y,2));
	double z = position.z;
	
	if(useThikDisk) {
		if(r<A1) {
			n += n1*cos(M_PI*r/(2*A1))/cos(M_PI*Rsun/(2*A1))/pow(cosh(z/H1),2);
		}
	}
	
	if(useThinDisk) {
		n += n2*exp(-pow(r-A2,2)/pow(A2,2))/pow(cosh(z/H2),2);
	}
	
	if(useSpiralArm) {
		double nA = 0;
		double gj;
		double sj;
		double xDist[1000];
		double yDist[1000];
		for(int j = 0; j<5;j++) {
			gj=0;
			for(int i = 0 ; i<1000;i++){
				xDist[i] = xArm[i][j] - x;
				yDist[i] = yArm[i][j] - y;
			}
			double minDistx = *std::min_element(xDist,xDist+1000);
			double minDisty = *std::min_element(yDist,yDist+1000);
			sj =sqrt( pow(minDistx,2)+pow(minDisty,2));
			
			if(r<Aa)
				gj = exp(-pow(sj/wjwa[j],2))/pow(cosh(r-Aa),2)/2;
			n += fjna[j]*gj/pow(cosh(z/hjHa[j]),2);
		}
	}
	
}

} //namespace
		
	

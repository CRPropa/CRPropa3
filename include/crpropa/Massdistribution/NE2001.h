#ifndef CRPROPA_NE2001_H
#define CRPROPA_NE2001_H

#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Massdistribution/Density.h"

#include <string>
#include <vector>
#include <math.h>

namespace crpropa {

class NE2001: public Density {

	bool isforHI = false;
	bool isforHII = true;
	bool isforH2 = false;
	
	bool useSpiralArm = true;
	bool useThikDisk = true;
	bool useThinDisk = true;
	bool useGalacticCenter = true;
//	bool useVoids = false;
//	bool useClumps = false;
//	bool useLISM = false;
	
	double Rsun;
	double n1, H1, A1;	//Parameter for Thik Disk
	double n2, H2, A2;	//Parameter for Thin Disk
	double fjna[5];		//SpiralArm Parameter
	double hjHa[5];
	double wjwa[5];
	double nGC, RGC, hGC; //Parameter for GalacticCenter
	double Aa;
	Vector3d PositionGC;
	double xArm[1000][5];	// x-Data for SpiralArmPosition
	double yArm[1000][5];	// y-Data for SpiralArmPosition
	
public:
	NE2001();	
	double getDensity(const Vector3d &position) const;
	double getHIIDensity(const Vector3d &position) const;
	
	bool getisforHI();
	bool getisforHII();
	bool getisforH2();
	
	bool getuseSpiralArm();
	bool getuseThikDisk();
	bool getuseThinDisk();
	bool getuseGalacticCenter();
/*	bool getuseVoids();
	bool getuseClumps();
	bool getuseLISM(); */
	
	void setuseSpiralArm(bool use);
	void setuseThikDisk(bool use);
	void setuseThinDisk(bool use);
	void setuseGalacticCenter(bool use);
/*	void setuseVoids(bool use);
	void setuseClumps(bool use);
	void setuseLISM(bool use); */
	
	
	
	
	
};

} //namespace crpropa

#endif //CRPROPA_NE2001_H

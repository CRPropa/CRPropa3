#ifndef CRPROPA_POHL2008_H
#define CRPROPA_POHL2008_H

#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Massdistribution/Density.h"
#include "crpropa/Grid.h"
#include "crpropa/GridTools.h"

#include <string>

namespace crpropa {

/**
	@class Pohl08 
	@brief Grid for H2 density based on arxiv:0712.4264 Fits avalibal at http://www.app.physik.uni-potsdam.de/gas.html 
*/

class Pohl08: public Density {
	ScalarGrid Grid;
		
	// Information on Modeltyp. DO NOT CHANGE
	bool isforHI = false;
	bool isforHII = false;
	bool isforH2 = true; 
	bool useReducedGrid=true;	// use every 4th value of Pohls Grid -> equal spacing in all axis 

public:
	Pohl08();
	void loadPohlGrid();
	double getDensity(const Vector3d &position) const;
	double getH2Density(const Vector3d &position) const;
	
	bool getisforHI();
	bool getisforHII();
	bool getisforH2();
	
	bool getuseReducedGrid();
	void setuseReducedGrid(bool reduced);
	
	
};


} //namespace crpropa

#endif //CRPROPA_POHL2008_H

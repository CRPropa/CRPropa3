#include "crpropa/Massdistribution/Pohl2008.h"

namespace crpropa {


void Pohl08::loadPohlGrid() {
	if(useReducedGrid)
	{
		Grid = ScalarGrid(Vector3d(-14950*pc,-14950*pc,-487.5*pc), 300,300,10,100*pc); 
		loadGrid(&Grid, "share/crpropa/Pohl_Reduced.txt", 2); //Conversion 2 for double Target in H2
	}
	else
	{
		Grid = ScalarGrid(Vector3d(-14950*pc,-14950*pc, -487.5*pc),300,300,40,100*pc);
		loadGrid(&Grid, "share/crpropa/Pohl.txt",2); //conversion 2 for double Target in H2
	}
}


double Pohl08::getDensity(const Vector3d &position) const {
	return getH2Density(position);
}

double Pohl08::getH2Density(const Vector3d &position) const {
	
	Vector3d pos = position;
	if(!useReducedGrid)
	{	//transformation of Position for spacing of grid
		pos.z *=4;
	}
	
	return Grid.interpolate(pos);
}

bool Pohl08::getisforHI() {
	return isforHI;
}

bool Pohl08::getisforHII() {
	return isforHII;
}

bool Pohl08::getisforH2() {
	return isforH2;
}

bool Pohl08::getuseReducedGrid() {
	return useReducedGrid;
}

void Pohl08::setuseReducedGrid(bool reduced) {

	useReducedGrid = reduced;
	loadPohlGrid();		//changing the dimension of Grid
}			

} //namespace

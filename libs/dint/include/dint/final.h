#ifndef _FINAL_H_
#define _FINAL_H_

#include "dint/spectrum.h"
#include "dint/cvector.h"

void CheckEnergy(const int sourceTypeSwitch, const double brightPhaseExp,
		 const double startingRedshift, 
		 const double rightRedshift, const Spectrum* pSpectrum, 
		 const dCVector* pEnergy, const double initialTotalEnergy);
void FinalPrintOutToTheScreen(const double distance, 
			      const double startingRedshift,
			      const double propagatingDistance);

#endif

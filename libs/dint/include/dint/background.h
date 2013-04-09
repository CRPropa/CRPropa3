
#ifndef _BACKGROUND_H_
#define _BACKGROUND_H_

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>

#include "dint/cvector.h"
#include "dint/error.h"
#include "dint/utilities.h"
#include "dint/const.h"
#include "dint/gauleg.h"

#define GAULEG_POINTS 31 // number of gaussian quadrature points for integration

using namespace std;

void LoadPhotonBackground(const double redshift, 
                          dCVector* pBgEnergy, dCVector* pBgEnergyWidth, 
                          dCVector* pBgPhotonDensity, const int aIRFlag,
			  const double aZmax_IR, const int aRadioFlag,
			  const double aH0, const double aOmegaM, const double aOmegaLambda);
void LoadCMB(const double redshift, const dCVector* pBgEnergy, 
             const dCVector* pBgEnergyWidth, dCVector* pBgPhotonDensity);
void LoadIR(const double redshift, const dCVector* pBgEnergy,
            const dCVector* pBgEnergyWidth, dCVector* pBgPhotonDensity,
	    const int aIRFlag, const double aZmax_IR);
double IR2(const double redshift, const double BgEnergy);
double HighIR(const double zTarget, const double zObserve, 
	      const double energy0, const double deltaO,
	      const double deltaD);
double LowIR(const double zTarget, const double zObserve, 
	     const double energy0, const double deltaO,
	     const double deltaD);
double OpticalIR(const double energy);
double DustIR(const double energy);
void LoadRadio(const double redshift, const dCVector* pBgEnergy,
               const dCVector* pBgEnergyWidth, dCVector* pBgPhotonDensity,
	       const int aRadioFlag, const double aH0, const double aOmegaM, 
	       const double aOmegaLambda);
double HighRadio(const double zTarget, const double zObserve, 
		 const double energy0, const double aH0, 
		 const double aOmegaM, const double aOmegaLambda);
double MedRadio(const double zTarget, const double zObserve, 
		const double energy0, const double aH0, const double aOmegaM, 
		const double aOmegaLambda);
double ObsRadio(const double zTarget, const double zObserve, 
		const double energy0,  const double aH0, 
		const double aOmegaM, const double aOmegaLambda);
/*
// Routine added by Gunter (July 2005) :
void DumpBgSpectrum(const dCVector* pBgEnergy, const dCVector* pBgEnergyWidth,
		    const dCVector* pBgPhotonDensity, const char* filename);
*/

#endif

#include <stdio.h>
#include <math.h>
#include "dint/spectrum.h"
#include "dint/const.h"

void CheckEnergy(const int sourceTypeSwitch, const double brightPhaseExp,
                 const double startingRedshift,
		 const double rightRedshift, const Spectrum* pSpectrum,
		 const dCVector* pEnergy, const double initialTotalEnergy)
{
    int i;
    double tempFraction;
    double photonNumber = 0.;
    double leptonNumber = 0.;
    double nucleonNumber = 0.;
    double neutrinoNumber = 0.;
    double totalNumber = 0.;
    double photonEnergy = 0.;
    double leptonEnergy = 0.;
    double nucleonEnergy = 0.;
    double neutrinoEnergy = 0.;
    double totalEnergy = 0.;
    double energyRedshiftFactor;

    tempFraction = (1. + startingRedshift)/(1. + rightRedshift);
    for (i = 0; i < pEnergy->dimension; i++)
    {
        photonNumber += (pSpectrum->spectrum)[PHOTON][i];
        photonEnergy += (pSpectrum->spectrum)[PHOTON][i]*(pEnergy->vector)[i];
        leptonNumber += (pSpectrum->spectrum)[ELECTRON][i] +
	    (pSpectrum->spectrum)[POSITRON][i];
        leptonEnergy += (pEnergy->vector)[i]*
	    ((pSpectrum->spectrum)[ELECTRON][i] +
            (pSpectrum->spectrum)[POSITRON][i]);
    }

    nucleonNumber = GetNucleonNumber(pSpectrum);
    nucleonEnergy = GetNucleonEnergy(pSpectrum, pEnergy);
    neutrinoNumber = GetNeutrinoNumber(pSpectrum);
    neutrinoEnergy = GetNeutrinoEnergy(pSpectrum, pEnergy);
    totalNumber = photonNumber + leptonNumber + nucleonNumber + neutrinoNumber;
    totalEnergy = photonEnergy + leptonEnergy + nucleonEnergy + neutrinoEnergy;

    if (sourceTypeSwitch == 0)    /* single source */
	energyRedshiftFactor = 1./pow(tempFraction, 4);
    else
    {
	energyRedshiftFactor = C_0/H_0*1.e6/1.e5*pow(1. + startingRedshift,
            -1.5)/pow(tempFraction,1.5+brightPhaseExp)/(brightPhaseExp-2.5)*
	    (pow(tempFraction,brightPhaseExp-2.5) - 1.);
    }

    //    printf("    Total energy = %15.6E (arbitrary units)\n", totalEnergy);
    //    printf("    (total energy)/(injected energy) = %g\n",
    //	   totalEnergy/initialTotalEnergy/energyRedshiftFactor);
}


void FinalPrintOutToTheScreen(const double distance,
			      const double startingRedshift,
			      const double propagatingDistance)
{
    double analyticalDistance;

    printf("\n\nTotal distance was %g Mpc ",
        distance/DISTANCE_UNIT/1.e6);
    analyticalDistance = 2./3.*C_0*1.e-5/H_0*(1. -
        pow(1. + startingRedshift, -3./2.));
    printf("vs. real distance %g Mpc.\n", analyticalDistance);
    printf("And total x was %g pc\n", propagatingDistance);
}

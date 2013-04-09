#include <math.h>
#include "dint/const.h"
#include "dint/utilities.h"


/* NOTE: note that these functions use MIN_ENERGY_EXP, etc. instead of
   NUC_MIN_ENERGY_EXP, etc.  Thus these functions are tied with the main
   routine only. */
double PionToPhoton(const int iPhoton, const int iPion)
{
    double piEnergyRight;
    double piEnergyLeft;
    double deltaPiE;
    double result;
    
    piEnergyRight = pow(10., MIN_ENERGY_EXP + (iPion+0.5)/BINS_PER_DECADE)/
        ELECTRON_MASS;
    piEnergyLeft = pow(10., MIN_ENERGY_EXP + (iPion-0.5)/BINS_PER_DECADE)/
        ELECTRON_MASS;
    deltaPiE = piEnergyRight - piEnergyLeft;
    
    if (iPhoton < iPion)
        result = 2./deltaPiE/BINS_PER_DECADE*log(10.);
    else if (iPhoton == iPion)
    {
        result = 2.*(1./deltaPiE - piEnergyLeft/deltaPiE/deltaPiE/
            BINS_PER_DECADE*log(10.));
    }
    else
        result = 0.;
        
    if (result < 0.)
	Error("PionToPhoton: decay spectrum negative.", PROGRAM_ERROR);
    
    return result;
}


double PionToLepton(const double leptonEnergy, const double pionEnergy)
{
    const double r = (1.056595e8/1.39567e8)*(1.056595e8/1.39567e8);
    const double A0 = 0.94486e0;
    const double A2 = -2.7892e0;
    const double A3 = 1.2397e0;
    const double B0 = -2.4126e0;
    const double B0p = -2.8951e0;
    const double B2 = 4.3426e0;
    const double B3 = -1.9300e0;
    
    double ratio;
    double result;

    ratio = leptonEnergy/pionEnergy;
    if (leptonEnergy < r*pionEnergy)
    {
        result = 1./(1.0 - r)/pionEnergy*(A0 + A2*ratio*ratio + A3*ratio*
            ratio*ratio);
    }
    else
    {
        result = 1./(1.0 - r)/pionEnergy*(B0 + B0p*log(ratio) + 
            B2*ratio*ratio + B3*ratio*ratio*ratio);
    }
    if (result < 0.) 
        result = 0.;
    
    return result;
}

double PionToElectronNeutrino(const double neutrinoEnergy, 
			      const double pionEnergy)
{
    const double r = (1.056595e8/1.39567e8)*(1.056595e8/1.39567e8);
    const double C0 = 1.1053e0;
    const double C2 = -4.46883e0;
    const double C3 = 3.71887e0;
    const double D0 = 13.846e0;
    const double D0p = 5.37053e0;
    const double D1 = -28.1116e0;
    const double D2 = 20.0558e0;
    const double D3 = -5.7902e0;

    double ratio;
    double result;

    ratio = neutrinoEnergy/pionEnergy;
    if (neutrinoEnergy < r*pionEnergy)
    {
	result = 1./(1. - r)/pionEnergy*(C0 + C2*ratio*ratio + C3*ratio*ratio
	    *ratio);
    }
    else
    {
	result = 1./(1. - r)/pionEnergy*(D0 + D0p*log(ratio) + D1*ratio +
            D2*ratio*ratio + D3*ratio*ratio*ratio);
    }
    if (result < 0.)
	result = 0.;

    return result;
}

double PionToMuonNeutrino(const int iNeutrino, const int iPion)
{
    const double r = (1.056595e8/1.39567e8)*(1.056595e8/1.39567e8);

    double neutEnergyLeft;
    double neutEnergyRight;
    double piEnergyLeft;
    double piEnergyRight;
    double deltaNeutE;
    double deltaPiE;
    double result;

    neutEnergyRight = pow(10., MIN_ENERGY_EXP + (iNeutrino+0.5)/
        BINS_PER_DECADE)/ELECTRON_MASS;
    neutEnergyLeft = pow(10., MIN_ENERGY_EXP + (iNeutrino-0.5)/
        BINS_PER_DECADE)/ELECTRON_MASS;
    deltaNeutE = neutEnergyRight - neutEnergyLeft;

    piEnergyRight = pow(10., MIN_ENERGY_EXP + (iPion+0.5)/BINS_PER_DECADE)/
        ELECTRON_MASS;
    piEnergyLeft = pow(10., MIN_ENERGY_EXP + (iPion-0.5)/BINS_PER_DECADE)/
        ELECTRON_MASS;
    deltaPiE = piEnergyRight - piEnergyLeft;

    if (neutEnergyRight < piEnergyLeft*(1.-r))
	result = log(piEnergyRight/piEnergyLeft)/(1. - r)/deltaPiE;
    else if (neutEnergyLeft < piEnergyLeft*(1.-r))
    {
	result = (deltaNeutE/(1. - r)*log(piEnergyRight*(1. - r)/
            neutEnergyRight) + neutEnergyRight/(1. - r) - piEnergyLeft -
            neutEnergyLeft/(1. - r)*log(neutEnergyRight/(piEnergyLeft*
            (1. - r))))/deltaPiE/deltaNeutE;
    }
    else if ((neutEnergyLeft < piEnergyRight*(1.-r)) && 
	     (neutEnergyRight > piEnergyRight*(1.-r)))
    {
	result = (piEnergyRight - neutEnergyLeft/(1. - r) - neutEnergyLeft/
            (1. - r)*log(piEnergyRight*(1. - r)/neutEnergyLeft))/deltaPiE/
            deltaNeutE;
    }
    else
	result = 0.;

    if (result < 0.)
	result = 0.;

    return result;
}

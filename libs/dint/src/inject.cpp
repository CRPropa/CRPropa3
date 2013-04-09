#include <math.h>
#include "dint/cvector.h"
#include "dint/spectrum.h"
#include "dint/const.h"
#include "dint/frag.h"
#include "dint/decay.h"
#include "dint/utilities.h"
#include "dint/inject.h"

// E.Armengaud - Dec 2005
// This routine is not used anymore :
//   the injection spectrum must now be computed within CRPropa.

void SetInjectionSpectrum(const PARTICLE part, const double InjEnergy, 
			  const double HInjEnergy, const double deltaE_hadron,
			  const dCVector* pEnergy,
			  const dCVector* pEnergyWidth, Spectrum* pQ_0) {
  int num_main_bins, maxBin, i;
  double criticalEnergy;
  
  InitializeSpectrum(pQ_0);
  num_main_bins = pEnergy->dimension;
  
  if ((pEnergyWidth->dimension != num_main_bins) ||
      (pQ_0->numberOfMainBins != num_main_bins))
    Error("PhotonMonoInjection: inconsistent dimensions", PROGRAM_ERROR);
  
  if (part != NOTHING) {
    criticalEnergy = InjEnergy/ELECTRON_MASS;
    maxBin = (int)((log10(criticalEnergy*ELECTRON_MASS) -
		    MAX_ENERGY_EXP)*BINS_PER_DECADE + num_main_bins);
    (pQ_0->spectrum)[part][maxBin] = 1.;
    
  } else {
    // In this case, we model the injection spectrum created by pair production
    // with a power law of index -7/4
    if (deltaE_hadron == 0.e0) Error("DeltaE_Hadron = 0 !", PROGRAM_ERROR);
    double sum=0.;
    criticalEnergy = HInjEnergy/ELECTRON_MASS;
    for (i = 0; i < num_main_bins; i++) {
      if (pEnergy->vector[i] < criticalEnergy) {
	(pQ_0->spectrum)[ELECTRON][i] = pow(pEnergy->vector[i],-7./4.)*
	  (pEnergyWidth->vector)[i];
	sum += (pQ_0->spectrum)[ELECTRON][i]*(pEnergy->vector)[i];
      }
    }
    sum *= ELECTRON_MASS;
    sum = deltaE_hadron/sum/2.;
    for (i = 0; i < num_main_bins; i++) {
      (pQ_0->spectrum)[ELECTRON][i] *= sum;
      (pQ_0->spectrum)[POSITRON][i] = (pQ_0->spectrum)[ELECTRON][i];
    }
    
  }
}

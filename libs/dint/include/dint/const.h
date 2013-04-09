
#ifndef _CONST_H_
#define _CONST_H_


#define ELECTRON_MASS 5.110e5
#define DISTANCE_UNIT 3.0856e18
#define VOLUME_UNIT 6.652448e-25*3.0856e18
#define PI 3.141592
#define C 3.e10
#define H_0 71.
// Added July 2005 : cosmological parameters
#define OMEGA_M 0.3// NOW THEY ARE NOT defined like this, but taken as input parameters.
#define OMEGA_LAMBDA 0.7 // if not specified as parameters, take these values

// CHANGE (Guenter; 7/20/1998)
#define DMAX  1.e6
#define CLUSTER_DISTANCE 100.
#define CLUSTER_FACTOR 1.
#define SOURCE_CLUSTER_DISTANCE 0.1
#define SOURCE_CLUSTER_FACTOR 1.

#define NUM_IP_ELEMENTS 1309125
#define NUM_IS_ELEMENTS 571200
#define NUM_PP_ELEMENTS 558050
#define NUM_TPP_ELEMENTS 548250
// photopion production
#define NUM_PPP_PROTON_SCAT_ELEMENTS 111541
#define NUM_PPP_PROTON_NEUTRON_ELEMENTS 111541
#define NUM_PPP_PROTON_PHOTON_ELEMENTS 20000
#define NUM_PPP_PROTON_ELECTRON_ELEMENTS 20000
#define NUM_PPP_PROTON_POSITRON_ELEMENTS 20000
// nucleon pair production
#define NUM_NPP_ELEMENTS 20000
// neutrino production from PPP
#define NUM_PPP_PROTON_ANTI_ELECTRON_NEUTRINO_ELEMENTS 20000
#define NUM_PPP_PROTON_MUON_NEUTRINO_ELEMENTS 20000
#define NUM_PPP_PROTON_ANTI_MUON_NEUTRINO_ELEMENTS 20000

// Parameters used in the code
// Warning : those cannot be changed unless various parts of CRPropa are also changed!
#define BINS_PER_DECADE 10    // number of bins per decade
#define MIN_ENERGY_EXP 7      // minimum spectrum energy = 10 MeV
#define MAX_ENERGY_EXP (MIN_ENERGY_EXP + 17) // maximum spectrum energy

#define BG_MIN_ENERGY_EXP (-8 - MIN_ENERGY_EXP + 7)
#define BG_MAX_ENERGY_EXP (2 - MIN_ENERGY_EXP + 7)
#define EM_MIN_ENERGY_EXP (MIN_ENERGY_EXP)
#define NUC_MIN_ENERGY_EXP (14 + MIN_ENERGY_EXP - 7)
#define NEUT_MIN_ENERGY_EXP (17 + MIN_ENERGY_EXP - 7)
#define NUM_MAIN_BINS ((MAX_ENERGY_EXP - MIN_ENERGY_EXP)*BINS_PER_DECADE)
#define NUM_BG_BINS ((BG_MAX_ENERGY_EXP - BG_MIN_ENERGY_EXP)*BINS_PER_DECADE)

#define EM_NUM_MAIN_BINS ((MAX_ENERGY_EXP - EM_MIN_ENERGY_EXP)*BINS_PER_DECADE)
#define NUC_NUM_MAIN_BINS ((MAX_ENERGY_EXP - NUC_MIN_ENERGY_EXP)*BINS_PER_DECADE)
#define NEUT_NUM_MAIN_BINS ((MAX_ENERGY_EXP - NEUT_MIN_ENERGY_EXP)*BINS_PER_DECADE)


#endif

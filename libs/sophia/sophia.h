#ifndef _SOPHIA_H_
#define _SOPHIA_H_

extern "C" {
void sophiaevent_(int&, double&, double[][2000], int[], int&, double&, int&,
		double&, int&, double[], double[]);
}

/*
 The arguments are the following (see sophiaevent.f) :

 - nature, Ein = input nature and energy of the nucleon
 nature = 0 -> p ; 1 -> n
 Ein : in GeV (Sophia standard energy unit)
 - OutPart,OutPartType,NbOutPart = output data, respectively :
 P(2000,5) list of 4-momenta + masses of output particles (in GeV)
 LList(2000) list of output particle IDs
 NP nb of output particles
 Particle IDs are :

 cc      13  proton
 cc      14  neutron
 cc      -13 antiproton
 cc      -14 antineutron
 cc      1   photon
 cc      2   e+
 cc      3   e-
 cc      15  nu_e
 cc      16  antinu_e
 cc      17  nu_muon
 cc      18  antinu_muon

 - z_particle : needed to estimate the CMB temperature or IRB density
 - bgFlag = 1 for CMB, 2 for Primack et al. (1999) IRB. NO OTHER REDSHIFT MODEL
 IS IMPLEMENTED, AND ITS REDSHIFT EVOLUTION IS "TRIVIAL", ie. same as CMB.
 - Zmax_IRB : the photon density of IRB is null above this redshift

 Warning with passing arguments (variables->pointers) in fortran->C !

 */

//#define sophiaevent sophiaevent_

#endif

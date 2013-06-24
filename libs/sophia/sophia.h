#ifndef _SOPHIA_H
#define _SOPHIA_H

extern "C" {
void sophiaevent_(int&, double&, double[][2000], int[], int&, double&, int&,
		double&, int&, double[], double[]);
}

/*
 The arguments are the following
 - nature: 0 -> p, 1 -> n
 - input energy of nucleon in GeV
 - list of 4-momenta + masses of output particles (in GeV)
	13  proton
	14  neutron
	-13 antiproton
	-14 antineutron
	1   photon
	2   e+
	3   e-
	15  nu_e
	16  antinu_e
	17  nu_muon
	18  antinu_muon
 - list of output particle ids
 - number of output particles
 - redshift
 - photon background flag: 1 -> CMB, 2 -> IRB Kneiske
 	 (Primack et al. (1999) IRB is outcommented in sophia_interface.f on line 16320
 - maximum redshift: the photon density of IRB is null above this redshift
 -
 -
 -
 */

#endif

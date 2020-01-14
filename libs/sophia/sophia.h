#ifndef _SOPHIA_H
#define _SOPHIA_H


extern "C" {
void sophiaevent_(int& nature,                   // IN:  channel: 0 -> p, 1 -> n
									double& Ein,                   // IN:  input energy of nucleon in GeV
									double& eps,                   // IN:  energy of target photon in GeV
									double outputEnergy[5][2000],  // OUT: list of 4-momenta + masses of output particles (in GeV)
									int outPartID[2000],           // OUT: list of output particle IDs (see list below)
									int& nParticles                // OUT: number of output particles
		);
}

/*
 - list of output particle IDs
		 13  proton
		-13  anti-proton
		 14  neutron
		-14  anti-neutron
		 1   photon
		 2   e+
		 3   e-
		 15  nu_e
		 16  anti-nu_e
		 17  nu_mu
		 18  anti-nu_mu
*/

#endif

#ifndef _SOPHIA_H
#define _SOPHIA_H

extern "C" {
/*
int type                    particle type: 0=proton, 1=neutron
double inputEnergy          input energy of nucleon in GeV
double eps                  energy of target photon in eV
double momentaList[2000]    list of output particle energies in GeV
double idList[2000]         list of output particle IDs
int nParticles              number of output particles

 - list of output particle ids
	 13  proton
	-13  antiproton
	 14  neutron
	-14  antineutron
	 1   photo
	 2   e+
	 3   e-
	 15  nu_e
	 16  antinu_e
	 17  nu_mu
	 18  antinu_mu
*/
void sophiaevent_(int& type, double& inputEnergy, double& eps,
				  double momentaList[2000], int idList[], int& nParticles);
}




#endif

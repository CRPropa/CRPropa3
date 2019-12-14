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
 - list of output particle ids
     13  proton
     14  neutron
    -13  antiproton
    -14  antineutron
     1   photon
     2   e+
     3   e-
     15  nu_e
     16  antinu_e
     17  nu_muon
     18  antinu_muon
*/

#endif

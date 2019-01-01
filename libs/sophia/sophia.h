#ifndef _SOPHIA_H
#define _SOPHIA_H


extern "C" {
void sophiaevent_(int& type, double& inputEnergy, double& e, double momentum[2000],
        int id[], int& n);
}


// void sophiaevent_(int& nature,  // channel: 0 -> p, 1 -> n
//                double& Ein,  // input energy of nucleon in GeV
//                double momentaList[][2000],  // list of 4-momenta + masses of output particles (in GeV)
//                double particleList[2000],
//                int& nParticles,  // number of output particles
//                double& eps  // energy of target photon in eV
//                );
// }


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

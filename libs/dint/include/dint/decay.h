#ifndef _DECAY_H
#define _DECAY_H

double PionToPhoton(const int iPhoton, const int iPion);
double PionToLepton(const double leptonEnergy, const double pionEnergy);
double PionToElectronNeutrino(const double neutrinoEnergy,
			      const double pionEnergy);
double PionToMuonNeutrino(const int iNeutrino, const int iPion);

#endif

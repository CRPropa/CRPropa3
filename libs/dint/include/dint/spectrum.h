#ifndef DINT__SPECTRUM_H
#define DINT__SPECTRUM_H

#include "dint/vector.h"
#include "dint/cvector.h"

#define NUM_SPECIES 11

typedef enum PARTICLE {
	PHOTON = 0,
	ELECTRON = 1,
	POSITRON = 2,
	PROTON = 3,
	NEUTRON = 4,
	ELECTRON_NEUTRINO = 5,
	ANTI_ELECTRON_NEUTRINO = 6,
	MUON_NEUTRINO = 7,
	ANTI_MUON_NEUTRINO = 8,
	TAU_NEUTRINO = 9,
	ANTI_TAU_NEUTRINO = 10,
	NOTHING = 11
} PARTICLE;
// NOTHING type added by Gunter (2005) for secondary from pair production

/* Spectrum is a struct that corresponds to the total spectrum.  It
   has all particle spectra in it.  Since particle data is more or less
   fixed and independent, we do not specify the number of particle species. */
typedef struct Spectrum {
	int numberOfMainBins;
	dMatrix spectrum;
} Spectrum;

void NewSpectrum(Spectrum* pSpectrum, const int num_bins);
void DeleteSpectrum(Spectrum* pSpectrum);
void InitializeSpectrum(Spectrum* pSpectrum);
void InitializeEMSpectrum(Spectrum* pSpectrum);
void InitializeNucleonSpectrum(Spectrum* pSpectrum);
void InitializeNeutrinoSpectrum(Spectrum* pSpectrum);
void AddSpectrum(Spectrum* a, const Spectrum* b);
void SetSpectrum(Spectrum* pSpectrum1, const Spectrum* pSpectrum2);
void SetEMSpectrum(Spectrum* pSpectrum1, const Spectrum* pSpectrum2);
void SetNucleonSpectrum(Spectrum* pSpectrum1, const Spectrum* pSpectrum2);
void SetNeutrinoSpectrum(Spectrum* pSpectrum1, const Spectrum* pSpectrum2);
double GetNumber(const Spectrum* pSpectrum);
double GetEMNumber(const Spectrum* pSpectrum);
double GetNucleonNumber(const Spectrum* pSpectrum);
double GetNeutrinoNumber(const Spectrum* pSpectrum);
double GetEnergy(const Spectrum* pSpectrum, const dCVector* pEnergy);
double GetEMEnergy(const Spectrum* pSpectrum, const dCVector* pEnergy);
double GetNucleonEnergy(const Spectrum* pSpectrum, const dCVector* pEnergy);
double GetNeutrinoEnergy(const Spectrum* pSpectrum, const dCVector* pEnergy);
void DumpSpectrum(const dCVector* pEnergy, const dCVector* pEnergyWidth,
		const Spectrum* pSpectrum, const char* filename);

#endif

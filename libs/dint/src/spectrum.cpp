#include <stdio.h>
#include <string.h>
#include "dint/spectrum.h"
#include "dint/vector.h"
#include "dint/cvector.h"
#include "dint/const.h"
#include "dint/utilities.h"


void NewSpectrum(Spectrum* pSpectrum, const int num_bins) {
  pSpectrum->numberOfMainBins = num_bins;
  pSpectrum->spectrum = New_dMatrix(NUM_SPECIES, num_bins);
  
  InitializeSpectrum(pSpectrum);    // always initialize!
}

void DeleteSpectrum(Spectrum* pSpectrum) {
  Delete_dMatrix(pSpectrum->spectrum);
}

void InitializeSpectrum(Spectrum* pSpectrum) {
  for (int i=0; i<NUM_SPECIES; i++) {
    for (int j=0; j<pSpectrum->numberOfMainBins; j++)
      pSpectrum->spectrum[i][j] = 0.;
  }
}

void InitializeEMSpectrum(Spectrum* pSpectrum) {
  for (int i=0; i<3; i++) {
    for (int j=0; j<pSpectrum->numberOfMainBins; j++)
      pSpectrum->spectrum[i][j] = 0.;
  }
}

void InitializeNucleonSpectrum(Spectrum* pSpectrum) {
  for (int i=3; i<5; i++) {
    for (int j=0; j<pSpectrum->numberOfMainBins; j++)
      pSpectrum->spectrum[i][j] = 0.;
  }
}

void InitializeNeutrinoSpectrum(Spectrum* pSpectrum) {
  for (int i=5; i<NUM_SPECIES; i++) {
    for (int j=0; j<pSpectrum->numberOfMainBins; j++)
      pSpectrum->spectrum[i][j] = 0.;
  }
}


void SetSpectrum(Spectrum* pSpectrum1, const Spectrum* pSpectrum2) {
  if (pSpectrum1->numberOfMainBins != pSpectrum2->numberOfMainBins)
    Error("SetSpectrum: inconsistent dimensions", PROGRAM_ERROR);
  
  for (int i=0; i<NUM_SPECIES; i++) {
    for (int j=0; j<pSpectrum1->numberOfMainBins; j++)
      pSpectrum1->spectrum[i][j] = pSpectrum2->spectrum[i][j];
  }
}

void SetEMSpectrum(Spectrum* pSpectrum1, const Spectrum* pSpectrum2) {
  if (pSpectrum1->numberOfMainBins != pSpectrum2->numberOfMainBins)
    Error("SetEMSpectrum: inconsistent dimension", PROGRAM_ERROR);

  for (int i=0; i<3; i++) {
    for (int j=0; j<pSpectrum1->numberOfMainBins; j++)
      pSpectrum1->spectrum[i][j] = pSpectrum2->spectrum[i][j];
  }
}

void SetNucleonSpectrum(Spectrum* pSpectrum1, const Spectrum* pSpectrum2) {
  if (pSpectrum1->numberOfMainBins != pSpectrum2->numberOfMainBins)
    Error("SetSpectrum: inconsistent dimension", PROGRAM_ERROR);
  
  for (int i=3; i<5; i++) {
    for (int j=0; j<pSpectrum1->numberOfMainBins; j++)
      pSpectrum1->spectrum[i][j] = pSpectrum2->spectrum[i][j];
  }
}

void SetNeutrinoSpectrum(Spectrum* pSpectrum1, const Spectrum* pSpectrum2) {
  if (pSpectrum1->numberOfMainBins != pSpectrum2->numberOfMainBins)
    Error("SetSpectrum: inconsistent dimension", PROGRAM_ERROR);

  for (int i=5; i<NUM_SPECIES; i++) {
    for (int j=0; j<pSpectrum1->numberOfMainBins; j++)
      pSpectrum1->spectrum[i][j] = pSpectrum2->spectrum[i][j];
  }
}

double GetNumber(const Spectrum* pSpectrum) {
  double number = 0.;
    
  for (int i=0; i<NUM_SPECIES; i++) {
    for (int j=0; j<pSpectrum->numberOfMainBins; j++)
      number += pSpectrum->spectrum[i][j];
  }
  return number;
}

double GetEMNumber(const Spectrum* pSpectrum) {
  double number = 0.;
    
  for (int i=0; i<3; i++) {
    for (int j=0; j<pSpectrum->numberOfMainBins; j++)
      number += pSpectrum->spectrum[i][j];
  }
  return number;
}

double GetNucleonNumber(const Spectrum* pSpectrum) {
  double number = 0.;
  
  for (int i=3; i<5; i++) {
    for (int j=0; j<pSpectrum->numberOfMainBins; j++)
      number += pSpectrum->spectrum[i][j];
  }
  return number;
}

double GetNeutrinoNumber(const Spectrum* pSpectrum) {
  double number = 0.;
  
  for (int i=5; i<NUM_SPECIES; i++) {
    for (int j=0; j<pSpectrum->numberOfMainBins; j++)
      number += pSpectrum->spectrum[i][j];
  }
  return number;
}

double GetEnergy(const Spectrum* pSpectrum, const dCVector* pEnergy) {
  double totalEnergy = 0.;
  
  if (pSpectrum->numberOfMainBins != pEnergy->dimension)
    Error("GetEnergy: inconsistent dimensions", PROGRAM_ERROR);
  
  for (int i=0; i<NUM_SPECIES; i++) {
    for (int j=0; j<pSpectrum->numberOfMainBins; j++)
      totalEnergy += (pSpectrum->spectrum)[i][j]*(pEnergy->vector)[j];
  }
  return totalEnergy;
}

double GetEMEnergy(const Spectrum* pSpectrum, const dCVector* pEnergy) {
  double totalEnergy = 0.;
  
  if (pSpectrum->numberOfMainBins != pEnergy->dimension)
    Error("GetEnergy: inconsistent dimensions", PROGRAM_ERROR);
  
  for (int i=0; i<3; i++) {
    for (int j=0; j<pSpectrum->numberOfMainBins; j++)
      totalEnergy += (pSpectrum->spectrum)[i][j]*(pEnergy->vector)[j];
  }
  return totalEnergy;
}

double GetNucleonEnergy(const Spectrum* pSpectrum, const dCVector* pEnergy) {
  double totalEnergy = 0.;

  if (pSpectrum->numberOfMainBins != pEnergy->dimension)
    Error("GetEnergy: inconsistent dimensions", PROGRAM_ERROR);

  for (int i=3; i<5; i++) {
    for (int j=0; j<pSpectrum->numberOfMainBins; j++)
      totalEnergy += (pSpectrum->spectrum)[i][j]*(pEnergy->vector)[j];
  }
  return totalEnergy;
}

double GetNeutrinoEnergy(const Spectrum* pSpectrum, const dCVector* pEnergy) {
  double totalEnergy = 0.;

  if (pSpectrum->numberOfMainBins != pEnergy->dimension)
    Error("GetEnergy: inconsistent dimensions", PROGRAM_ERROR);

  for (int i=5; i<NUM_SPECIES; i++) {
    for (int j=0; j<pSpectrum->numberOfMainBins; j++)
      totalEnergy += (pSpectrum->spectrum)[i][j]*(pEnergy->vector)[j];
  }
  return totalEnergy;
}

void DumpSpectrum(const dCVector* pEnergy, const dCVector* pEnergyWidth,
		  const Spectrum* pSpectrum, const char* filename) {
  FILE* dumpFile;
  char f1[80] = "datafiles/";
  
  if ((pEnergy->dimension != pEnergyWidth->dimension) ||
      (pEnergyWidth->dimension != pSpectrum->numberOfMainBins))
    Error("DumpSpectrum: inconsistent dimensions", PROGRAM_ERROR);
  
  strncat(f1, filename, 79 - strlen(filename));
  dumpFile = SafeFOpen(f1, "w");
  // this is to send the dump file to a different directory (datafiles)
  //   by Guenter (7/20/1998)
  
  for (int i=0; i<NUM_SPECIES; i++) {
    for (int j=0; j<pSpectrum->numberOfMainBins; j++) {
      fprintf(dumpFile, "%15.4E %15.4E\n", 
	      ELECTRON_MASS*(pEnergy->vector)[j], 
	      (pSpectrum->spectrum)[i][j]/(pEnergyWidth->vector)[j]*
	      (pEnergy->vector)[j]*(pEnergy->vector)[j]);
      // proper unit conversion for energy
    }
  }
  fclose(dumpFile);
}

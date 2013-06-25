#ifndef DINT__LOAD_H
#define DINT__LOAD_H

#include "dint/rate.h"
#include "dint/utilities.h"
#include "dint/const.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "binfread.h"

using namespace std;

void LoadICSTables(RawTotalRate* ICSTotalRate, RawDiffRate* ICSPhotonRate,
                   RawDiffRate* ICSScatRate, const int num_main_bins,
		   string aDirTables);
void LoadPPTables(RawTotalRate* PPTotalRate, RawDiffRate* PPDiffRate, 
		  const int num_main_bins,
		  string aDirTables);
void LoadTPPTables(RawTotalRate* TPPTotalRate, RawDiffRate* TPPDiffRate, 
		   const int num_main_bins, string aDirTables); 
void LoadDPPTables(RawTotalRate* DPPRate, const int num_main_bins,
		   string aDirTables);
void LoadPPPNucleonTables(RawTotalRate* PPPProtonLossRate, 
			  RawTotalRate* PPPNeutronLossRate,
			  RawDiffRate* PPPProtonScatRate,
			  RawDiffRate* PPPProtonNeutronRate,
			  RawDiffRate* PPPNeutronProtonRate,
			  const int num_main_bins, string aDirTables);
void LoadPPPEMTables(RawDiffRate* PPPProtonPhotonRate,
		     RawDiffRate* PPPProtonElectronRate,
		     RawDiffRate* PPPProtonPositronRate,
		     RawDiffRate* PPPNeutronElectronRate, 
		     const int num_main_bins, string aDirTables);
void LoadNPPNucleonTables(RawTotalRate* NPPTotalRate, const int num_main_bins, 
			  string aDirTables);
void LoadNPPSecondaryTables(RawDiffRate* NPPDiffRate, const int num_main_bins, 
			    string aDirTables);
void LoadPPPNeutrinoTables(RawDiffRate* PPPProtonElectronNeutrinoRate,
			   RawDiffRate* PPPProtonAntiElectronNeutrinoRate,
			   RawDiffRate* PPPProtonMuonNeutrinoRate,
			   RawDiffRate* PPPProtonAntiMuonNeutrinoRate,
			   RawDiffRate* PPPNeutronAntiElectronNeutrinoRate,
			   RawDiffRate* PPPNeutronMuonNeutrinoRate,
			   RawDiffRate* PPPNeutronAntiMuonNeutrinoRate, 
			   const int num_main_bins, string aDirTables);
void LoadNeutronDecayNucleonTables(TotalRate* neutronDecayRate, 
				   DiffRate* neutronDecayProtonRate,
				   const int num_main_bins, string aDirTables);
void LoadNeutronDecaySecondaryTables(DiffRate* neutronDecayElectronRate,
				     const int num_main_bins, 
				     string aDirTables);
void LoadNeutrinoTables(const int tauNeutrinoMassSwitch,
			TotalRate* NNElNeutTotalRate, 
			TotalRate* NNMuonNeutTotalRate,
			TotalRate* NNTauNeutTotalRate, 
			DiffRate* NNElNeutScatRate, 
			DiffRate* NNElNeutMuonNeutRate,
			DiffRate* NNElNeutTauNeutRate, 
			DiffRate* NNElNeutElectronRate, 
			DiffRate* NNElNeutPhotonRate,
			DiffRate* NNElNeutProtonRate, 
			DiffRate* NNMuonNeutScatRate, 
			DiffRate* NNMuonNeutElNeutRate,
			DiffRate* NNMuonNeutTauNeutRate, 
			DiffRate* NNMuonNeutElectronRate, 
			DiffRate* NNMuonNeutPhotonRate, 
			DiffRate* NNMuonNeutProtonRate, 
			DiffRate* NNTauNeutScatRate,
			DiffRate* NNTauNeutElNeutRate, 
			DiffRate* NNTauNeutMuonNeutRate, 
			DiffRate* NNTauNeutElectronRate, 
			DiffRate* NNTauNeutPhotonRate, 
			DiffRate* NNTauNeutProtonRate, 
			const int num_main_bins, string aDirTables);

#endif

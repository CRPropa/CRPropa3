
#include "dint/load.h"

// Modified Apr 2005 : aDirTable argument

void LoadICSTables(RawTotalRate* ICSTotalRate, RawDiffRate* ICSPhotonRate,
		   RawDiffRate* ICSScatRate, const int num_main_bins,
		   string aDirTables)
{
  ReadRawTotalRate(ICSTotalRate, (aDirTables+"/ICSLoss.dat").c_str());
  ModifyRawTotalRate(ICSTotalRate, num_main_bins);

  ReadRawDiffRate(ICSPhotonRate, (aDirTables+"/ICSLP.dat").c_str());
  ModifyRawDiffRate(ICSPhotonRate, num_main_bins);

  ReadRawDiffRate(ICSScatRate, (aDirTables+"/ICSLS.dat").c_str());
  ModifyRawDiffRate(ICSScatRate, num_main_bins);
}


void LoadPPTables(RawTotalRate* PPTotalRate, RawDiffRate* PPDiffRate, 
		  const int num_main_bins, string aDirTables)
{
  ReadRawTotalRate(PPTotalRate, (aDirTables+"/PPLoss.dat").c_str());
  ModifyRawTotalRate(PPTotalRate, num_main_bins);
  
  ReadRawDiffRate(PPDiffRate, (aDirTables+"/PPPL.dat").c_str());
  ModifyRawDiffRate(PPDiffRate, num_main_bins);
}


void LoadTPPTables(RawTotalRate* TPPTotalRate, RawDiffRate* TPPDiffRate, 
		   const int num_main_bins, string aDirTables)
{
  ReadRawTotalRate(TPPTotalRate, (aDirTables+"/TPPLoss.dat").c_str());
  ModifyRawTotalRate(TPPTotalRate, num_main_bins);
  
  ReadRawDiffRate(TPPDiffRate, (aDirTables+"/TPPDiff.dat").c_str());
  ModifyRawDiffRate(TPPDiffRate, num_main_bins);
}

void LoadDPPTables(RawTotalRate* DPPRate, const int num_main_bins,
		   string aDirTables)
{
  ReadRawTotalRate(DPPRate, (aDirTables+"/DPP.dat").c_str());
  ModifyRawTotalRate(DPPRate, num_main_bins);
}

void LoadPPPNucleonTables(RawTotalRate* PPPProtonLossRate, 
			  RawTotalRate* PPPNeutronLossRate,
			  RawDiffRate* PPPProtonScatRate,
			  RawDiffRate* PPPProtonNeutronRate,
			  RawDiffRate* PPPNeutronProtonRate,
			  const int num_main_bins,
			  string aDirTables)
{
  FILE* PPPLoss;
  
  PPPLoss = SafeFOpen((aDirTables+"/PPPLoss.dat").c_str(), "r");
  binfread(PPPProtonLossRate->totalRate[0], sizeof(double), 
	(PPPProtonLossRate->mainDimension)*(PPPProtonLossRate->bgDimension), 
        PPPLoss);
  binfread(PPPNeutronLossRate->totalRate[0], sizeof(double), 
	(PPPNeutronLossRate->mainDimension)*(PPPNeutronLossRate->bgDimension), 
        PPPLoss);
  fclose(PPPLoss);
  ModifyRawTotalRate(PPPProtonLossRate, num_main_bins);
  ModifyRawTotalRate(PPPNeutronLossRate, num_main_bins);
  /* this is treated a little differently because the file contains both
     tables */


  ReadRawDiffRate(PPPProtonScatRate, (aDirTables+"/PPPPrS.dat").c_str());
  ModifyRawDiffRate(PPPProtonScatRate, num_main_bins);
  
  ReadRawDiffRate(PPPProtonNeutronRate, (aDirTables+"/PPPPrN.dat").c_str());
  ModifyRawDiffRate(PPPProtonNeutronRate, num_main_bins);
  
  ReadRawDiffRate(PPPNeutronProtonRate, (aDirTables+"/PPPNPr.dat").c_str());
  ModifyRawDiffRate(PPPNeutronProtonRate, num_main_bins);
}


void LoadPPPEMTables(RawDiffRate* PPPProtonPhotonRate,
		     RawDiffRate* PPPProtonElectronRate,
		     RawDiffRate* PPPProtonPositronRate,
		     RawDiffRate* PPPNeutronElectronRate, 
		     const int num_main_bins,
		     string aDirTables)
{
  FILE* PPPNE;
  
  ReadRawDiffRate(PPPProtonPhotonRate, (aDirTables+"/PPPPrPh.dat").c_str());
  ModifyRawDiffRate(PPPProtonPhotonRate, num_main_bins);
  
  ReadRawDiffRate(PPPProtonElectronRate, (aDirTables+"/PPPPrE.dat").c_str());
  ModifyRawDiffRate(PPPProtonElectronRate, num_main_bins);
  
  ReadRawDiffRate(PPPProtonPositronRate, (aDirTables+"/PPPPrPos.dat").c_str());
  
  PPPNE = SafeFOpen((aDirTables+"/PPPNE.dat").c_str(), "r");
  if (PPPNeutronElectronRate->numberOfElements != 
      PPPProtonPositronRate->numberOfElements)
    {
      Error("LoadPPPTables: inconsistent dimensions", PROGRAM_ERROR);
    }
  binfread(PPPNeutronElectronRate->diffRate, sizeof(double), 
	PPPProtonPositronRate->numberOfElements, PPPNE);
  fclose(PPPNE);
  CopyRawDiffRate(PPPNeutronElectronRate, PPPProtonPositronRate);
  ModifyRawDiffRate(PPPProtonPositronRate, num_main_bins);
  ModifyRawDiffRate(PPPNeutronElectronRate, num_main_bins);
}


void LoadNPPNucleonTables(RawTotalRate* NPPTotalRate, const int num_main_bins,
			  string aDirTables)
{
  ReadRawTotalRate(NPPTotalRate, (aDirTables+"/NPPLoss.dat").c_str());
  ModifyRawTotalRate(NPPTotalRate, num_main_bins);
}

void LoadNPPSecondaryTables(RawDiffRate* NPPDiffRate, const int num_main_bins,
			    string aDirTables)
{
  ReadRawDiffRate(NPPDiffRate, (aDirTables+"/NPPDiff.dat").c_str());
  ModifyRawDiffRate(NPPDiffRate, num_main_bins);
}

void LoadPPPNeutrinoTables(RawDiffRate* PPPProtonElectronNeutrinoRate,
			   RawDiffRate* PPPProtonAntiElectronNeutrinoRate,
			   RawDiffRate* PPPProtonMuonNeutrinoRate,
			   RawDiffRate* PPPProtonAntiMuonNeutrinoRate,
			   RawDiffRate* PPPNeutronAntiElectronNeutrinoRate,
			   RawDiffRate* PPPNeutronMuonNeutrinoRate,
			   RawDiffRate* PPPNeutronAntiMuonNeutrinoRate, 
			   const int num_main_bins, string aDirTables)
{
  FILE* PPPPrEN;
  FILE* PPPNAEN;
  FILE* PPPNMN;
  FILE* PPPNAMN;
  

  /* these 4 rates share the same bound */
  ReadRawDiffRate(PPPProtonAntiMuonNeutrinoRate, (aDirTables+"/PPPPrAMN.dat").c_str());
  
  if (PPPProtonElectronNeutrinoRate->numberOfElements != 
      PPPProtonAntiMuonNeutrinoRate->numberOfElements)
    {
      Error("LoadPPPTables: inconsistent dimensions", PROGRAM_ERROR);
    }
  PPPPrEN = SafeFOpen((aDirTables+"/PPPPrEN.dat").c_str(), "r");
  binfread(PPPProtonElectronNeutrinoRate->diffRate, sizeof(double),
	PPPProtonAntiMuonNeutrinoRate->numberOfElements, PPPPrEN);
  fclose(PPPPrEN);
  CopyRawDiffRateBound(PPPProtonElectronNeutrinoRate,
		       PPPProtonAntiMuonNeutrinoRate);
  
  if (PPPNeutronAntiElectronNeutrinoRate->numberOfElements != 
      PPPProtonAntiMuonNeutrinoRate->numberOfElements)
    {
      Error("LoadPPPTables: inconsistent dimensions", PROGRAM_ERROR);
    }
  PPPNAEN = SafeFOpen((aDirTables+"/PPPNAEN.dat").c_str(), "r");
  binfread(PPPNeutronAntiElectronNeutrinoRate->diffRate, sizeof(double),
	PPPProtonAntiMuonNeutrinoRate->numberOfElements, PPPNAEN);
  fclose(PPPNAEN);
  CopyRawDiffRateBound(PPPNeutronAntiElectronNeutrinoRate,
		       PPPProtonAntiMuonNeutrinoRate);
  
  if (PPPNeutronMuonNeutrinoRate->numberOfElements != 
      PPPProtonAntiMuonNeutrinoRate->numberOfElements)
    {
      Error("LoadPPPTables: inconsistent dimensions", PROGRAM_ERROR);
    }
  PPPNMN = SafeFOpen((aDirTables+"/PPPNMN.dat").c_str(), "r");
  binfread(PPPNeutronMuonNeutrinoRate->diffRate, sizeof(double),
	PPPProtonAntiMuonNeutrinoRate->numberOfElements, PPPNMN);
  fclose(PPPNMN);
  CopyRawDiffRateBound(PPPNeutronMuonNeutrinoRate,
		       PPPProtonAntiMuonNeutrinoRate);
  
  ModifyRawDiffRate(PPPProtonAntiMuonNeutrinoRate, num_main_bins);
  ModifyRawDiffRate(PPPProtonElectronNeutrinoRate, num_main_bins);
  ModifyRawDiffRate(PPPNeutronAntiElectronNeutrinoRate, num_main_bins);
  ModifyRawDiffRate(PPPNeutronMuonNeutrinoRate, num_main_bins);
  
  
  
  ReadRawDiffRate(PPPProtonAntiElectronNeutrinoRate, 
		  (aDirTables+"/PPPPrAEN.dat").c_str());
  ModifyRawDiffRate(PPPProtonAntiElectronNeutrinoRate, num_main_bins);
  
  /* these 2 share bound */
  ReadRawDiffRate(PPPProtonMuonNeutrinoRate, (aDirTables+"/PPPPrMN.dat").c_str());
  
  if (PPPNeutronAntiMuonNeutrinoRate->numberOfElements != 
      PPPProtonMuonNeutrinoRate->numberOfElements)
    {
      Error("LoadPPPTables: inconsistent dimensions", PROGRAM_ERROR);
    }
  PPPNAMN = SafeFOpen((aDirTables+"/PPPNAMN.dat").c_str(), "r");
  binfread(PPPNeutronAntiMuonNeutrinoRate->diffRate, sizeof(double),
	PPPProtonMuonNeutrinoRate->numberOfElements, PPPNAMN);
  fclose(PPPNAMN);
  CopyRawDiffRateBound(PPPNeutronAntiMuonNeutrinoRate,
		       PPPProtonMuonNeutrinoRate);
  
  ModifyRawDiffRate(PPPProtonMuonNeutrinoRate, num_main_bins);
  ModifyRawDiffRate(PPPNeutronAntiMuonNeutrinoRate, num_main_bins);
}

void LoadNeutronDecayNucleonTables(TotalRate* neutronDecayRate,
				   DiffRate* neutronDecayProtonRate,
				   const int num_main_bins,
				   string aDirTables)
{
  ReadTotalRate(neutronDecayRate, (aDirTables+"/neutronDecayLoss.dat").c_str());
  ModifyTotalRate(neutronDecayRate, num_main_bins);
  
  ReadDiffRate(neutronDecayProtonRate, (aDirTables+"/neutronDecayProton.dat").c_str());
  ModifyDiffRate(neutronDecayProtonRate, num_main_bins);
}

void LoadNeutronDecaySecondaryTables(DiffRate* neutronDecayElectronRate, 
				     const int num_main_bins,
				     string aDirTables)
{
  ReadDiffRate(neutronDecayElectronRate, (aDirTables+"/neutronDecayElectron.dat").c_str());
  ModifyDiffRate(neutronDecayElectronRate, num_main_bins);
}

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
			const int num_main_bins,
			string aDirTables)
{
  if (tauNeutrinoMassSwitch == 0)
    {
      ReadTotalRate(NNElNeutTotalRate, (aDirTables+"/ENT_0.dat").c_str());
      ReadTotalRate(NNMuonNeutTotalRate, (aDirTables+"/MNT_0.dat").c_str());
      ReadTotalRate(NNTauNeutTotalRate, (aDirTables+"/TNT_0.dat").c_str());
      
      ReadDiffRate(NNElNeutScatRate, (aDirTables+"/ENS_0.dat").c_str());
      ReadDiffRate(NNElNeutMuonNeutRate, (aDirTables+"/ENMN_0.dat").c_str());
      ReadDiffRate(NNElNeutTauNeutRate, (aDirTables+"/ENTN_0.dat").c_str());
      ReadDiffRate(NNElNeutElectronRate, (aDirTables+"/ENE_0.dat").c_str());
      ReadDiffRate(NNElNeutPhotonRate, (aDirTables+"/ENPh_0.dat").c_str());
      ReadDiffRate(NNElNeutProtonRate, (aDirTables+"/ENPr_0.dat").c_str());
      ReadDiffRate(NNMuonNeutElNeutRate, (aDirTables+"/MNEN_0.dat").c_str());
      ReadDiffRate(NNMuonNeutScatRate, (aDirTables+"/MNS_0.dat").c_str());
      ReadDiffRate(NNMuonNeutTauNeutRate, (aDirTables+"/MNTN_0.dat").c_str());
      ReadDiffRate(NNMuonNeutElectronRate, (aDirTables+"/MNE_0.dat").c_str());
      ReadDiffRate(NNMuonNeutPhotonRate, (aDirTables+"/MNPh_0.dat").c_str());
      ReadDiffRate(NNMuonNeutProtonRate, (aDirTables+"/MNPr_0.dat").c_str());
      ReadDiffRate(NNTauNeutElNeutRate, (aDirTables+"/TNEN_0.dat").c_str());
      ReadDiffRate(NNTauNeutMuonNeutRate, (aDirTables+"/TNMN_0.dat").c_str());
      ReadDiffRate(NNTauNeutScatRate, (aDirTables+"/TNS_0.dat").c_str());
      ReadDiffRate(NNTauNeutElectronRate, (aDirTables+"/TNE_0.dat").c_str());
      ReadDiffRate(NNTauNeutPhotonRate, (aDirTables+"/TNPh_0.dat").c_str());
      ReadDiffRate(NNTauNeutProtonRate, (aDirTables+"/TNPr_0.dat").c_str());
    }
  else
    {
      ReadTotalRate(NNElNeutTotalRate, (aDirTables+"/ENT_10.dat").c_str());
      ReadTotalRate(NNMuonNeutTotalRate, (aDirTables+"/MNT_10.dat").c_str());
      ReadTotalRate(NNTauNeutTotalRate, (aDirTables+"/TNT_10.dat").c_str());
      ReadDiffRate(NNElNeutScatRate, (aDirTables+"/ENS_10.dat").c_str());
      ReadDiffRate(NNElNeutMuonNeutRate, (aDirTables+"/ENMN_10.dat").c_str());
      ReadDiffRate(NNElNeutTauNeutRate, (aDirTables+"/ENTN_10.dat").c_str());
      ReadDiffRate(NNElNeutElectronRate, (aDirTables+"/ENE_10.dat").c_str());
      ReadDiffRate(NNElNeutPhotonRate, (aDirTables+"/ENPh_10.dat").c_str());
      ReadDiffRate(NNElNeutProtonRate, (aDirTables+"/ENPr_10.dat").c_str());
      ReadDiffRate(NNMuonNeutElNeutRate, (aDirTables+"/MNEN_10.dat").c_str());
      ReadDiffRate(NNMuonNeutScatRate, (aDirTables+"/MNS_10.dat").c_str());
      ReadDiffRate(NNMuonNeutTauNeutRate, (aDirTables+"/MNTN_10.dat").c_str());
      ReadDiffRate(NNMuonNeutElectronRate, (aDirTables+"/MNE_10.dat").c_str());
      ReadDiffRate(NNMuonNeutPhotonRate, (aDirTables+"/MNPh_10.dat").c_str());
      ReadDiffRate(NNMuonNeutProtonRate, (aDirTables+"/MNPr_10.dat").c_str());
      ReadDiffRate(NNTauNeutElNeutRate, (aDirTables+"/TNEN_10.dat").c_str());
      ReadDiffRate(NNTauNeutMuonNeutRate, (aDirTables+"/TNMN_10.dat").c_str());
      ReadDiffRate(NNTauNeutScatRate, (aDirTables+"/TNS_10.dat").c_str());
      ReadDiffRate(NNTauNeutElectronRate, (aDirTables+"/TNE_10.dat").c_str());
      ReadDiffRate(NNTauNeutPhotonRate, (aDirTables+"/TNPh_10.dat").c_str());
      ReadDiffRate(NNTauNeutProtonRate, (aDirTables+"/TNPr_10.dat").c_str());
    }

  ModifyTotalRate(NNElNeutTotalRate, num_main_bins);
  ModifyTotalRate(NNMuonNeutTotalRate, num_main_bins);
  ModifyTotalRate(NNTauNeutTotalRate, num_main_bins);
  ModifyDiffRate(NNElNeutScatRate, num_main_bins);
  ModifyDiffRate(NNElNeutMuonNeutRate, num_main_bins);
  ModifyDiffRate(NNElNeutTauNeutRate, num_main_bins);
  ModifyDiffRate(NNElNeutElectronRate, num_main_bins);
  ModifyDiffRate(NNElNeutPhotonRate, num_main_bins);
  ModifyDiffRate(NNElNeutProtonRate, num_main_bins);
  ModifyDiffRate(NNMuonNeutElNeutRate, num_main_bins);
  ModifyDiffRate(NNMuonNeutScatRate, num_main_bins);
  ModifyDiffRate(NNMuonNeutTauNeutRate, num_main_bins);
  ModifyDiffRate(NNMuonNeutElectronRate, num_main_bins);
  ModifyDiffRate(NNMuonNeutPhotonRate, num_main_bins);
  ModifyDiffRate(NNMuonNeutProtonRate, num_main_bins);
  ModifyDiffRate(NNTauNeutElNeutRate, num_main_bins);
  ModifyDiffRate(NNTauNeutMuonNeutRate, num_main_bins);
  ModifyDiffRate(NNTauNeutScatRate, num_main_bins);
  ModifyDiffRate(NNTauNeutElectronRate, num_main_bins);
  ModifyDiffRate(NNTauNeutPhotonRate, num_main_bins);
  ModifyDiffRate(NNTauNeutProtonRate, num_main_bins);
}

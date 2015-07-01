#ifndef DINT_PROPAGATION_H
#define DINT_PROPAGATION_H

#include "dint/rate.h"
#include "dint/const.h"
#include "dint/spectrum.h"
#include "dint/cvector.h"
#include "dint/load.h"
#include "dint/prepare.h"
#include "dint/sync.h"
#include "dint/inject.h"
#include "dint/background.h"
#include "dint/fold.h"
#include "dint/advance.h"
#include "dint/final.h"
#include "dint/utilities.h"


// DintEMCascade.
// Class based on the original DINT prop_second function, intended as starting
// point to create a simplified EM cascade calculation based on transport
// equations to replace dint completly in a future release of CRPropa.
class DintEMCascade
{
	private:
		//-------- Declaration of main variables --------
	  //---- Interaction table coefficients ----
	  RawTotalRate ICSTotalRate;
	  RawTotalRate PPTotalRate;
	  RawTotalRate TPPTotalRate;
	  RawTotalRate DPPRate;
	  
	  RawTotalRate PPPProtonLossRate;
	  RawTotalRate PPPNeutronLossRate;
	  RawTotalRate NPPTotalRate;
	  // total (interaction) rates before being folded into the background 
	
	  RawDiffRate ICSPhotonRate;
	  RawDiffRate ICSScatRate;
	  RawDiffRate PPDiffRate;
	  RawDiffRate TPPDiffRate;
	  
	  RawDiffRate PPPProtonScatRate;
	  RawDiffRate PPPProtonNeutronRate;
	  RawDiffRate PPPNeutronProtonRate;
	  RawDiffRate PPPProtonPhotonRate;
	  RawDiffRate PPPProtonElectronRate;
	  RawDiffRate PPPProtonPositronRate;
	  RawDiffRate PPPNeutronElectronRate;
	  RawDiffRate NPPDiffRate;
	  RawDiffRate PPPProtonElectronNeutrinoRate;
	  RawDiffRate PPPProtonAntiElectronNeutrinoRate;
	  RawDiffRate PPPProtonMuonNeutrinoRate;
	  RawDiffRate PPPProtonAntiMuonNeutrinoRate;
	  RawDiffRate PPPNeutronAntiElectronNeutrinoRate;
	  RawDiffRate PPPNeutronMuonNeutrinoRate;
	  RawDiffRate PPPNeutronAntiMuonNeutrinoRate;
	  // differential rates before being folded into the background
	
	  TotalRate neutronDecayRate;
	  DiffRate neutronDecayElectronRate;
	  DiffRate neutronDecayProtonRate;
	  
	  TotalRate NNElNeutTotalRate;
	  TotalRate NNMuonNeutTotalRate;
	  TotalRate NNTauNeutTotalRate;
	  // These total rates are net rates; i.e. scattered flux into same bin is subtracted
	  DiffRate NNElNeutScatRate;
	  DiffRate NNElNeutMuonNeutRate;
	  DiffRate NNElNeutTauNeutRate;
	  DiffRate NNElNeutElectronRate;
	  DiffRate NNElNeutPhotonRate;
	  DiffRate NNElNeutProtonRate;
	  DiffRate NNMuonNeutScatRate;
	  DiffRate NNMuonNeutElNeutRate;
	  DiffRate NNMuonNeutTauNeutRate;
	  DiffRate NNMuonNeutElectronRate;
	  DiffRate NNMuonNeutPhotonRate;
	  DiffRate NNMuonNeutProtonRate;
	  DiffRate NNTauNeutScatRate;
	  DiffRate NNTauNeutElNeutRate;
	  DiffRate NNTauNeutMuonNeutRate;
	  DiffRate NNTauNeutElectronRate;
	  DiffRate NNTauNeutPhotonRate;
	  DiffRate NNTauNeutProtonRate;
	  // rates from neutrino-neutrino interaction
		DiffRate syncRate;

		// Energy Bins
		dCVector deltaG; // dg used in continuous energy loss calculation

		dCVector bgEnergy;
		dCVector bgEnergyWidth;
		dCVector bgPhotonDensity;

		Spectrum Q_0;    // standard injection function
		Spectrum spectrumNew;
		Spectrum derivative;

	 //---- interaction rates folded with photon background ----
		TotalRate leptonTotalRate;
		TotalRate photonTotalRate;
		TotalRate protonTotalRate;
		TotalRate neutronTotalRate;
		
		DiffRate leptonScatRate;
		DiffRate leptonExchRate;
		DiffRate leptonPhotonRate;
		DiffRate photonLeptonRate;
		
		DiffRate protonScatRate;
		DiffRate protonNeutronRate;
		DiffRate neutronProtonRate;
		DiffRate protonPhotonRate;
		DiffRate protonElectronRate;
		DiffRate protonPositronRate;
		DiffRate neutronElectronRate;
		DiffRate neutronPositronRate;
		DiffRate protonElectronNeutrinoRate;
		DiffRate protonAntiElectronNeutrinoRate;
		DiffRate protonMuonNeutrinoRate;
		DiffRate protonAntiMuonNeutrinoRate;
		DiffRate neutronAntiElectronNeutrinoRate;
		DiffRate neutronMuonNeutrinoRate;
		DiffRate neutronAntiMuonNeutrinoRate;
		
		TotalRate elNeutTotalRate;
		TotalRate muonNeutTotalRate;
		TotalRate tauNeutTotalRate;
		
		DiffRate elNeutElectronRate;
		DiffRate elNeutPhotonRate;
		DiffRate muonNeutElectronRate;
		DiffRate muonNeutPhotonRate;
		DiffRate tauNeutElectronRate;
		DiffRate tauNeutPhotonRate;
		// rates from neutrino-neutrino interaction 
		
		dCVector synchrotronLoss;    // sgdot 
		dCVector otherLoss;    // tgdot 
		dCVector continuousLoss;   // gdot 

		dCVector pEnergy;
		dCVector pEnergyWidth;

		dCVector RedshiftArray ;
		dCVector DistanceArray ;

		// switches
		const int synchrotronSwitch;
		const int sourceTypeSwitch;  
		const int tauNeutrinoMassSwitch; 
		const int ICSSwitch;
		const int PPSwitch;
		const int TPPSwitch;
		const int DPPSwitch;
		const int PPPSwitch;
		const int NPPSwitch;
		const int neutronDecaySwitch;
		const int nucleonToSecondarySwitch;
		const int neutrinoNeutrinoSwitch;
		
		const int aIRFlag;
		const int aRadioFlag;

		const double aZmax_IR;

		const double aH0;
		const double aOmegaM;
		const double aOmegaLambda;

		dCVector pB_field;
		string aDirTables;
	public:

	DintEMCascade(int _aIRFlag, int _aRadioFlag, string _aDirTables, double B =
			1E-9, double _aH0 = H_0, double _aOmegaM = OMEGA_M, double _aOmegaLambda
			= OMEGA_LAMBDA);

	~DintEMCascade();
	
	void propagate(const double start_distance,
			const double stop_distance,
		 Spectrum* apInjectionSpectrum,
		 Spectrum* pSpectrum,
		 const double aCutcascade_Magfield = 0);
	 
};


#endif //DINT_PROPAGATION_H



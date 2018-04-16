#include "dint/prop_second.h"
#include "dint/DintEMCascade.h"

#include <string.h>
#include <stdlib.h>
#include <iostream>


DintEMCascade::DintEMCascade(int _aIRFlag, int _aRadioFlag, string _aDirTables,
		double B, double _aH0, double _aOmegaM, double
		_aOmegaLambda) : synchrotronSwitch(1), sourceTypeSwitch(0),
	tauNeutrinoMassSwitch(1), ICSSwitch(1), PPSwitch(1), TPPSwitch(1),
	DPPSwitch(1), PPPSwitch(0), NPPSwitch(0), neutronDecaySwitch(0),
	nucleonToSecondarySwitch(0), neutrinoNeutrinoSwitch(0), aIRFlag(_aIRFlag),
	aRadioFlag(_aRadioFlag), aH0(_aH0), aOmegaM(_aOmegaM),
	aOmegaLambda(_aOmegaLambda), aZmax_IR(5.), aDirTables(_aDirTables)
{
	BuildRedshiftTable(aH0, aOmegaM, aOmegaLambda, &RedshiftArray, &DistanceArray) ;

	// Initialize the bField
	New_dCVector(&pB_field, 1);
	for (size_t i = 0; i < 1; i++)
		pB_field.vector[i] = B; // B in gauss

	//-------- Set up energy bins --------
	New_dCVector(&deltaG, NUM_MAIN_BINS);
	New_dCVector(&bgEnergy, NUM_BG_BINS);
	New_dCVector(&bgEnergyWidth, NUM_BG_BINS);
	New_dCVector(&bgPhotonDensity, NUM_BG_BINS);

	New_dCVector(&pEnergy, NUM_MAIN_BINS);
	New_dCVector(&pEnergyWidth, NUM_MAIN_BINS);
	SetEnergyBins(MIN_ENERGY_EXP, &pEnergy, &pEnergyWidth);

	// set energy bins
	SetDeltaG(&pEnergy, &deltaG);

	SetEnergyBins(BG_MIN_ENERGY_EXP, &bgEnergy, &bgEnergyWidth);

	NewSpectrum(&Q_0, NUM_MAIN_BINS);
	NewSpectrum(&spectrumNew, NUM_MAIN_BINS);
	NewSpectrum(&derivative, NUM_MAIN_BINS);

	New_dCVector(&synchrotronLoss, NUM_MAIN_BINS);
	New_dCVector(&otherLoss, NUM_MAIN_BINS);
	New_dCVector(&continuousLoss, NUM_MAIN_BINS);

	if (ICSSwitch == 1) {
		NewRawTotalRate(&ICSTotalRate, EM_NUM_MAIN_BINS, NUM_BG_BINS);
		NewRawDiffRate(&ICSPhotonRate, EM_NUM_MAIN_BINS, NUM_BG_BINS,
			 NUM_IP_ELEMENTS);
		NewRawDiffRate(&ICSScatRate, EM_NUM_MAIN_BINS, NUM_BG_BINS,
			 NUM_IS_ELEMENTS);
	}

	if (PPSwitch == 1) {
		NewRawTotalRate(&PPTotalRate, EM_NUM_MAIN_BINS, NUM_BG_BINS);
		NewRawDiffRate(&PPDiffRate, EM_NUM_MAIN_BINS, NUM_BG_BINS,
			 NUM_PP_ELEMENTS);
	}

	if (TPPSwitch == 1) {
		NewRawTotalRate(&TPPTotalRate, EM_NUM_MAIN_BINS, NUM_BG_BINS);
		NewRawDiffRate(&TPPDiffRate, EM_NUM_MAIN_BINS, NUM_BG_BINS,
			 NUM_TPP_ELEMENTS);
	}

	if (DPPSwitch == 1)
		NewRawTotalRate(&DPPRate, EM_NUM_MAIN_BINS, NUM_BG_BINS);

	//---- read in coefficient tables; clipping is done here if necessary ----
	if (ICSSwitch == 1)
		LoadICSTables(&ICSTotalRate, &ICSPhotonRate, &ICSScatRate,
			NUM_MAIN_BINS, aDirTables);
	if (PPSwitch == 1)
		LoadPPTables(&PPTotalRate, &PPDiffRate, NUM_MAIN_BINS, aDirTables);
	if (TPPSwitch == 1)
		LoadTPPTables(&TPPTotalRate, &TPPDiffRate, NUM_MAIN_BINS, aDirTables);
	if (DPPSwitch == 1)
		LoadDPPTables(&DPPRate, NUM_MAIN_BINS, aDirTables);

	NewTotalRate(&leptonTotalRate, NUM_MAIN_BINS);
	NewTotalRate(&photonTotalRate, NUM_MAIN_BINS);

	NewTotalRate(&protonTotalRate, NUM_MAIN_BINS);
	NewTotalRate(&neutronTotalRate, NUM_MAIN_BINS);

	NewDiffRate(&leptonScatRate, NUM_MAIN_BINS);
	NewDiffRate(&leptonExchRate, NUM_MAIN_BINS);
	NewDiffRate(&leptonPhotonRate, NUM_MAIN_BINS);
	NewDiffRate(&photonLeptonRate, NUM_MAIN_BINS);

	NewDiffRate(&protonScatRate, NUM_MAIN_BINS);
	NewDiffRate(&protonNeutronRate, NUM_MAIN_BINS);
	NewDiffRate(&neutronProtonRate, NUM_MAIN_BINS);
	NewDiffRate(&protonPhotonRate, NUM_MAIN_BINS);
	NewDiffRate(&protonElectronRate, NUM_MAIN_BINS);
	NewDiffRate(&protonPositronRate, NUM_MAIN_BINS);
	NewDiffRate(&neutronElectronRate, NUM_MAIN_BINS);
	NewDiffRate(&neutronPositronRate, NUM_MAIN_BINS);
	NewDiffRate(&protonElectronNeutrinoRate, NUM_MAIN_BINS);
	NewDiffRate(&protonAntiElectronNeutrinoRate, NUM_MAIN_BINS);
	NewDiffRate(&protonMuonNeutrinoRate, NUM_MAIN_BINS);
	NewDiffRate(&protonAntiMuonNeutrinoRate, NUM_MAIN_BINS);
	NewDiffRate(&neutronAntiElectronNeutrinoRate, NUM_MAIN_BINS);
	NewDiffRate(&neutronMuonNeutrinoRate, NUM_MAIN_BINS);
	NewDiffRate(&neutronAntiMuonNeutrinoRate, NUM_MAIN_BINS);

	NewTotalRate(&elNeutTotalRate, NUM_MAIN_BINS);
	NewTotalRate(&muonNeutTotalRate, NUM_MAIN_BINS);
	NewTotalRate(&tauNeutTotalRate, NUM_MAIN_BINS);

	NewDiffRate(&elNeutElectronRate, NUM_MAIN_BINS);
	NewDiffRate(&elNeutPhotonRate, NUM_MAIN_BINS);
	NewDiffRate(&muonNeutElectronRate, NUM_MAIN_BINS);
	NewDiffRate(&muonNeutPhotonRate, NUM_MAIN_BINS);
	NewDiffRate(&tauNeutElectronRate, NUM_MAIN_BINS);
	NewDiffRate(&tauNeutPhotonRate, NUM_MAIN_BINS);
}

DintEMCascade::~DintEMCascade()
{
	DeleteDiffRate(&leptonScatRate);
	DeleteDiffRate(&leptonExchRate);
	DeleteDiffRate(&leptonPhotonRate);
	DeleteDiffRate(&photonLeptonRate);
	DeleteDiffRate(&protonScatRate);
	DeleteDiffRate(&protonNeutronRate);
	DeleteDiffRate(&neutronProtonRate);
	DeleteDiffRate(&protonPhotonRate);
	DeleteDiffRate(&protonElectronRate);
	DeleteDiffRate(&protonPositronRate);
	DeleteDiffRate(&neutronElectronRate);
	DeleteDiffRate(&neutronPositronRate);
	DeleteDiffRate(&protonElectronNeutrinoRate);
	DeleteDiffRate(&protonAntiElectronNeutrinoRate);
	DeleteDiffRate(&protonMuonNeutrinoRate);
	DeleteDiffRate(&protonAntiMuonNeutrinoRate);
	DeleteDiffRate(&neutronAntiElectronNeutrinoRate);
	DeleteDiffRate(&neutronMuonNeutrinoRate);
	DeleteDiffRate(&neutronAntiMuonNeutrinoRate);

	DeleteDiffRate(&elNeutElectronRate);
	DeleteDiffRate(&elNeutPhotonRate);
	DeleteDiffRate(&muonNeutElectronRate);
	DeleteDiffRate(&muonNeutPhotonRate);
	DeleteDiffRate(&tauNeutElectronRate);
	DeleteDiffRate(&tauNeutPhotonRate);

	DeleteTotalRate(&leptonTotalRate);
	DeleteTotalRate(&photonTotalRate);
	DeleteTotalRate(&protonTotalRate);
	DeleteTotalRate(&neutronTotalRate);
	DeleteTotalRate(&elNeutTotalRate);
	DeleteTotalRate(&muonNeutTotalRate);
	DeleteTotalRate(&tauNeutTotalRate);

	if (ICSSwitch == 1) {
		DeleteRawDiffRate(&ICSPhotonRate);
		DeleteRawDiffRate(&ICSScatRate);
		DeleteRawTotalRate(&ICSTotalRate);
	}
	if (PPSwitch == 1) {
		DeleteRawDiffRate(&PPDiffRate);
		DeleteRawTotalRate(&PPTotalRate);
	}
	if (TPPSwitch == 1) {
		DeleteRawDiffRate(&TPPDiffRate);
		DeleteRawTotalRate(&TPPTotalRate);
	}
	if (DPPSwitch == 1)
		DeleteRawTotalRate(&DPPRate);

	DeleteSpectrum(&Q_0);
	DeleteSpectrum(&spectrumNew);
	DeleteSpectrum(&derivative);

	Delete_dCVector(&synchrotronLoss);
	Delete_dCVector(&continuousLoss);
	Delete_dCVector(&otherLoss);

	Delete_dCVector(&deltaG);
	Delete_dCVector(&bgEnergy);
	Delete_dCVector(&bgEnergyWidth);
	Delete_dCVector(&bgPhotonDensity);

	Delete_dCVector(&RedshiftArray) ;
	Delete_dCVector(&DistanceArray) ;

	Delete_dCVector(&pEnergy);
	Delete_dCVector(&pEnergyWidth);
}


void DintEMCascade::propagate(const double start_distance,
		const double stop_distance, Spectrum* apInjectionSpectrum,
		Spectrum* pSpectrum, const double aCutcascade_Magfield) {

	double convergeParameter = 1.e-8;

	// -------- Redshift Estimation ---------------
	double startingRedshift = getRedshift(RedshiftArray, DistanceArray, start_distance) ;
	double stopRedshift = getRedshift(RedshiftArray, DistanceArray, stop_distance) ;

	double minDistance = 0.;
	double brightPhaseExp = 0.;

	//---- Initialize distance ----

	SetSpectrum(&Q_0, apInjectionSpectrum) ;
	PrepareSpectra(sourceTypeSwitch, &Q_0, pSpectrum, &spectrumNew, &derivative);

	//--------- START of actual computation --------
	//---- initialize indices and parameters ----
	double leftRedshift = startingRedshift;
	double deltaRedshift = 1. - (pEnergy.vector)[2]/(pEnergy.vector)[3];
	int lastIndex = 0;
	double propagatingDistance = 0.;
	double distance = 0.;  // distance variable: dist in FORTRAN

	//-------- This is where propagation takes place --------
	do {
		double rightRedshift, centralRedshift;
		// this firstIndex & lastIndex pair is used to bin redshift finer
		//  near the end of propagation (z close to 0)
		ComputeRedshifts(sourceTypeSwitch, leftRedshift, &deltaRedshift,
				 &rightRedshift, &centralRedshift, &lastIndex);
		// limit the distance step to the stop redshift
		if (rightRedshift < stopRedshift)
			rightRedshift = stopRedshift;

		//---- compute various distance parameters ----
		double redshiftRatio = (1. + rightRedshift)/(1. + leftRedshift);

		// (cosmological parameters added July 2005)
		double leftDistance = getDistance(RedshiftArray,DistanceArray,leftRedshift) ;
		double rightDistance =  getDistance(RedshiftArray,DistanceArray,rightRedshift) ;
		double distanceStep = leftDistance - rightDistance ;
		distance += distanceStep;

		int numSmallSteps = (int)(ceil(distanceStep/DISTANCE_UNIT/DMAX));
		double smallDistanceStep = distanceStep/DISTANCE_UNIT/numSmallSteps;
		double x = 0.;

		//---- compute the photon background at given redshift ----
		LoadPhotonBackground(centralRedshift, &bgEnergy, &bgEnergyWidth,
			 &bgPhotonDensity, aIRFlag, aZmax_IR, aRadioFlag,
			 aH0, aOmegaM, aOmegaLambda);

		//---- initialize rates ----
		InitializeLeptonCoefficients(&leptonTotalRate, &leptonScatRate,
				 &leptonExchRate, &leptonPhotonRate);
		InitializePhotonCoefficients(&photonTotalRate, &photonLeptonRate);

		Initialize_dCVector(&continuousLoss);

		Initialize_dCVector(&otherLoss);

		//---- fold interaction rates w/ photon background ----
		if (ICSSwitch == 1)
			FoldICS(&bgPhotonDensity, &ICSTotalRate, &ICSPhotonRate,
				&ICSScatRate, &leptonTotalRate, &leptonPhotonRate,
				&leptonScatRate);
		if (TPPSwitch == 1)
			FoldTPP(&bgPhotonDensity, &pEnergy, &TPPTotalRate, &TPPDiffRate,
				&leptonTotalRate, &leptonScatRate, &leptonExchRate,
				&otherLoss);
		if (PPSwitch == 1)
			FoldPP(&bgPhotonDensity, &PPTotalRate, &PPDiffRate,
			 &photonTotalRate, &photonLeptonRate);
		if (DPPSwitch == 1)
			FoldDPP(&bgPhotonDensity, &DPPRate, &photonTotalRate,
				&photonLeptonRate);

		//---- main iteration (convergence) block ----

		for (int loopCounter = 0; loopCounter < numSmallSteps; loopCounter++)
		{
			double B_loc = 0;
			if (synchrotronSwitch == 1)
			{
				NewDiffRate(&syncRate, NUM_MAIN_BINS);
				int B_bins = pB_field.dimension;
				int B_bin = (int)((double)(B_bins)*(propagatingDistance+x)/1.e6/
						start_distance);
				B_loc=(pB_field.vector)[B_bin];
				InitializeSynchrotron(B_loc, &pEnergy, &pEnergyWidth,
						&synchrotronLoss, &syncRate,
						aDirTables);
			}
			//---- compute continuous energy loss for electrons ----
			ComputeContinuousEnergyLoss(synchrotronSwitch, &synchrotronLoss,
					&otherLoss, &continuousLoss);

			// used for source evolution, disabled here
			double bkgFactor = 1.;
			double evolutionFactor = 0;

			AdvanceEMStep(sourceTypeSwitch, PPSwitch, ICSSwitch,
					TPPSwitch, DPPSwitch, synchrotronSwitch, PPPSwitch,
					NPPSwitch, neutronDecaySwitch, nucleonToSecondarySwitch,
					neutrinoNeutrinoSwitch, smallDistanceStep, evolutionFactor,
					convergeParameter, bkgFactor, &Q_0, &photonLeptonRate,
					&protonElectronRate, &neutronPositronRate,
					&protonPositronRate, &neutronElectronRate,
					&neutronDecayElectronRate, &elNeutElectronRate,
					&muonNeutElectronRate, &tauNeutElectronRate,
					&protonPhotonRate, &elNeutPhotonRate, &muonNeutPhotonRate,
					&tauNeutPhotonRate, &leptonTotalRate, &leptonScatRate,
					&leptonExchRate, &continuousLoss, &deltaG, &photonTotalRate,
					&leptonPhotonRate, &syncRate, pSpectrum, &spectrumNew);

			SetEMSpectrum(pSpectrum, &spectrumNew);
			// update spectrum

			if (aCutcascade_Magfield != 0 && B_loc != 0 )
			{
				// Estimate the effect of B field on the 1D approximation (added E.A. June 2006)
				bool lEcFlag = 0 ;
				int lIndex =  0;

				double a_ics = (3.-log10(4.))/4. ;
				double b_ics = pow(10.,8.-7.*a_ics) ;
				while (!lEcFlag)
				{
					double lEnergy = (pEnergy.vector)[lIndex] ;
					// Time scales are computed in parsec
					double t_sync = 3.84e6/(lEnergy*B_loc*B_loc*ELECTRON_MASS) ;
					double t_larmor = (1.1e-21)*ELECTRON_MASS*lEnergy/(B_loc*aCutcascade_Magfield) ;
					double t_ics;
					if (lEnergy <= 1.e15/ELECTRON_MASS) {
						t_ics = 4.e2*1.e15/(ELECTRON_MASS*lEnergy) ;
					} else if (lEnergy <= 1.e18/ELECTRON_MASS) {
						t_ics = 4.e2*lEnergy*ELECTRON_MASS/1.e15 ;
					} else if (lEnergy <= 1.e22/ELECTRON_MASS) {
						t_ics = b_ics*pow(lEnergy*ELECTRON_MASS/1.e15,a_ics) ;
					} else t_ics = 1.e8*lEnergy*ELECTRON_MASS/1.e22 ;
					if (t_larmor >= t_sync || t_larmor >= t_ics) lEcFlag = 1 ;
					// defines the "critical" energy : the e+/- spectrum is set to 0 for E<E_c
					(pSpectrum->spectrum)[ELECTRON][lIndex]=0 ;
					(pSpectrum->spectrum)[POSITRON][lIndex]=0 ;
					lIndex++ ;
				}
			}

			if (synchrotronSwitch == 1)     // synchrotron == true (on)
				DeleteDiffRate(&syncRate);

			x += smallDistanceStep;
		} // for loopCounter < numSmallSteps

		propagatingDistance += x;   // increment distance

		//---- redshift bins down ----
		// Force redshiftdown to use more accurate method, see issue  #174
		RedshiftDown(-1, redshiftRatio, &pEnergy, pSpectrum, &spectrumNew);

		//---- prepare for new step ----
		leftRedshift = rightRedshift;
	} while ((lastIndex != 1) && (leftRedshift > stopRedshift));

}


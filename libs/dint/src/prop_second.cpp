
#include "dint/prop_second.h"

#ifdef DEBUG
void DumpEnergy(const Spectrum* pSpectrum, const dCVector* pEnergy) {
  printf("Initial energy: %15.6E\n", GetEnergy(pSpectrum, pEnergy));
  printf("Nucleon fraction: %15.6E\n", 
	 GetNucleonEnergy(pSpectrum, pEnergy)/GetEnergy(pSpectrum, pEnergy));
  printf("electromagnetic fraction: %15.6E\n", 
	 GetEMEnergy(pSpectrum, pEnergy)/GetEnergy(pSpectrum, pEnergy));
}
#endif

//--------------------------------------------------------------------------
// Prop_second : see the file prop_second.h for explanations
//--------------------------------------------------------------------------

// Switches to read 1) the folded interaction rates : total pair production and 
//   total energy losses by photons; 2) the photon background density
//#define EXTRACT_PAIR_PROD_RATE
//#define EXTRACT_PHOTON_TOTAL_RATE
//#define PRINT_PHOTON_BACKGROUND
//#define EXTRACT_LEPTON_RATE

// Parameters of the redshift table (used for distance-redshift conversion)
#define NBINS_REDSHIFTTABLE 1000
#define ZMIN_REDSHIFTTABLE 0.0001
#define ZMAX_REDSHIFTTABLE 100

void prop_second(const double dist_observer,
		 //const double InjEnergy,
		 //const PARTICLE part, 
		 //const double HInjEnergy, const double deltaE_hadron,
		 const dCVector* pB_field,
		 const dCVector* pEnergy,
		 const dCVector* pEnergyWidth,
		 Spectrum* apInjectionSpectrum,
		 Spectrum* pSpectrum,
		 string aDirTables,
		 const int aIRFlag, const double aZmax_IR, const int aRadioFlag,
		 const double aH0 = H_0, const double aOmegaM = OMEGA_M, 
		 const double aOmegaLambda = OMEGA_LAMBDA,
		 const double aCutcascade_Magfield = 0) {
clock_t start = clock();
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


  //---- Main parameter variables ----
  int B_bin, B_bins;     // bin numbers for B-field array
  
  dCVector deltaG; // dg used in continuous energy loss calculation

  dCVector bgEnergy;
  dCVector bgEnergyWidth;
  dCVector bgPhotonDensity;
  
  double distance;    // main distance variable (cm) 

  int numSmallSteps;  // number of small steps in one redshift step
  double convergeParameter;   // redshift increment

  int synchrotronSwitch;
  double B_loc;     // local strength of extragalactic magnetic field (gauss)
  DiffRate syncRate;
  
  int sourceTypeSwitch;   // source type: single (0) or diffuse (1)
  double brightPhaseExp;  // bright phase exponent 
  double startingRedshift;    // zsource in FORTRAN 
  double minDistance;     // minimal distance to sources 
  double bkgFactor;       // local background enhancement factor 
  double totalEnergyInput;    // einj0 in FORTRAN 

  Spectrum Q_0;    // standard injection function
  Spectrum spectrumNew;
  Spectrum derivative;

  double initialPhotonEnergy; // totpe in FORTRAN 
  double initialLeptonEnergy; // totee 
  double initialNucleonEnergy;
  double initialNeutrinoEnergy;
  double initialTotalEnergy;  // tote 
  double initialPhotonNumber; // totpn 
  double initialLeptonNumber; // toten 
  double initialNucleonNumber;
  double initialNeutrinoNumber;
  double initialTotalNumber;  // totn 
  
  int iterationCounter;   // niter in FORTRAN 
  int firstIndex; // ifirst 
  double leftRedshift;    // zi 
  double deltaRedshift;   //zdelt 
  int lastIndex;  // ilast 
  double propagatingDistance; 
  // totalx: total propagating distance (pc) 
  double rightRedshift;   // zf 
  double centralRedshift; // zt 
  double redshiftRatio;   // erat 
  //  double tempCoefficient;     // coeff 
  double distanceStep;    // tzdist*tdz (cm) 
  double rightDistance;   // d1 (Mpc)
  double leftDistance;    // d2 (Mpc) 
  double evolutionFactor,evolutionFactor1; // zfact 
  double smallDistanceStep;   // dx (pc) 
  double x;   // x (pc) 

  // Energy below which the 1D approximation is not a priori correct :
  bool lEcFlag ;
  int lIndex ;
  double lEnergy, t_sync, t_larmor, t_ics = 0 ;
  double a_ics = (3.-log10(4.))/4. ;
  double b_ics = pow(10.,8.-7.*a_ics) ;
  // (used if the flag aCutcascade_Magfield is set)
  
  //    FILE* input;
  int tauNeutrinoMassSwitch;

  // interaction switches
  int ICSSwitch;
  int PPSwitch;
  int TPPSwitch;
  int DPPSwitch;
  
  int PPPSwitch;
  int NPPSwitch;
  int neutronDecaySwitch;
  int nucleonToSecondarySwitch;
  
  int neutrinoNeutrinoSwitch;
  
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
  
  DiffRate elNeutScatRate;
  DiffRate elNeutMuonNeutRate;
  DiffRate elNeutTauNeutRate;
  DiffRate elNeutElectronRate;
  DiffRate elNeutPhotonRate;
  DiffRate elNeutProtonRate;
  DiffRate muonNeutScatRate;
  DiffRate muonNeutElNeutRate;
  DiffRate muonNeutTauNeutRate;
  DiffRate muonNeutElectronRate;
  DiffRate muonNeutPhotonRate;
  DiffRate muonNeutProtonRate;
  DiffRate tauNeutScatRate;
  DiffRate tauNeutElNeutRate;
  DiffRate tauNeutMuonNeutRate;
  DiffRate tauNeutElectronRate;
  DiffRate tauNeutPhotonRate;
  DiffRate tauNeutProtonRate;
  // rates from neutrino-neutrino interaction 
  
  dCVector synchrotronLoss;    // sgdot 
  dCVector otherLoss;    // tgdot 
  dCVector continuousLoss;   // gdot 
  dCVector protonContinuousLoss;    // pgdot 

  int loopCounter;
  //-------- End of variable declaration --------

  //    numSmallSteps = 100;
  convergeParameter = 1.e-8;

  // -------- Redshift Estimation ---------------
  //  startingRedshift = pow(1. - 3.*H_0*1.e5*dist_observer/2./C, -2./3.) - 1.;
  // This was the old analytic redshift computation.
  dCVector RedshiftArray ;
  dCVector DistanceArray ;
  BuildRedshiftTable(aH0, aOmegaM, aOmegaLambda, &RedshiftArray, &DistanceArray) ;
  startingRedshift = getRedshift(RedshiftArray, DistanceArray, dist_observer) ;
  //    printf("distance/Mpc: %15.6E\n", dist_observer);
  //    printf("particle type: %d\n", part);
  //    printf("injection energy/eV: %15.6E\n", InjEnergy);

  B_bins = pB_field->dimension;
  
  synchrotronSwitch = 1;
  tauNeutrinoMassSwitch = 1;
  ICSSwitch = 1;
  PPSwitch = 1;
  TPPSwitch = 1;
  DPPSwitch = 1;
  //    PPPSwitch = 1;
  PPPSwitch = 0;
  NPPSwitch = 0;
  //     synchrotronSwitch = 0;
  //      tauNeutrinoMassSwitch = 0;
  //      ICSSwitch = 0;
  //      PPSwitch = 0;
  //      TPPSwitch = 0;
  //     DPPSwitch = 0;
      //    PPPSwitch = 1;
      PPPSwitch = 0;
      NPPSwitch = 0;
      //      printf("DINT modified!!\n");
#ifdef EXTRACT_PAIR_PROD_RATE
  NPPSwitch = 1; // (allows to compute new pair prod rate)
#endif

  //    neutronDecaySwitch = 1;
  neutronDecaySwitch = 0;
  nucleonToSecondarySwitch = 0;
  neutrinoNeutrinoSwitch = 1;
  sourceTypeSwitch = 0;
  minDistance = 0.;
  brightPhaseExp = 0.;
  
  //-------- Set up energy bins --------
  New_dCVector(&deltaG, NUM_MAIN_BINS);
  New_dCVector(&bgEnergy, NUM_BG_BINS);
  New_dCVector(&bgEnergyWidth, NUM_BG_BINS);
  New_dCVector(&bgPhotonDensity, NUM_BG_BINS);
  
  // set energy bins 
  SetDeltaG(pEnergy, &deltaG);
  
  SetEnergyBins(BG_MIN_ENERGY_EXP, &bgEnergy, &bgEnergyWidth);
  
  NewSpectrum(&Q_0, NUM_MAIN_BINS);
  NewSpectrum(&spectrumNew, NUM_MAIN_BINS);
  NewSpectrum(&derivative, NUM_MAIN_BINS);
  
  New_dCVector(&synchrotronLoss, NUM_MAIN_BINS);
  New_dCVector(&otherLoss, NUM_MAIN_BINS);
  New_dCVector(&continuousLoss, NUM_MAIN_BINS);
  New_dCVector(&protonContinuousLoss, NUM_MAIN_BINS);
    
    
  //---- Select injection model ----
  //  SetInjectionSpectrum(part, InjEnergy, HInjEnergy, deltaE_hadron,
  //		       pEnergy, pEnergyWidth, &Q_0);
  //    SetInjectionSpectrum(part, pEnergy, pEnergyWidth, &Q_0, InjEnergy);

  SetSpectrum(&Q_0, apInjectionSpectrum) ;
  // No call anymore to SetInjectionSpectrum, which is not useful anymore.
  totalEnergyInput = GetEnergy(&Q_0, pEnergy);
#ifdef DEBUG
  DumpEnergy(&Q_0, pEnergy);
#endif
  
  //-------- Create arrays --------
  // NOTE: I first make them table size; they will be "clipped" after reading in tables
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

  if (PPPSwitch == 1) {
    NewRawTotalRate(&PPPProtonLossRate, NUC_NUM_MAIN_BINS, NUM_BG_BINS);
    NewRawTotalRate(&PPPNeutronLossRate, NUC_NUM_MAIN_BINS, NUM_BG_BINS);
    
    NewRawDiffRate(&PPPProtonScatRate, NUC_NUM_MAIN_BINS, NUM_BG_BINS, 
		   NUM_PPP_PROTON_SCAT_ELEMENTS);
    NewRawDiffRate(&PPPProtonNeutronRate, NUC_NUM_MAIN_BINS, NUM_BG_BINS, 
		   NUM_PPP_PROTON_NEUTRON_ELEMENTS);
    NewRawDiffRate(&PPPNeutronProtonRate, NUC_NUM_MAIN_BINS, NUM_BG_BINS, 
		   NUM_PPP_PROTON_NEUTRON_ELEMENTS);
    
    if (nucleonToSecondarySwitch == 1) {
      NewRawDiffRate(&PPPProtonPhotonRate, NUC_NUM_MAIN_BINS, 
		     NUM_BG_BINS, NUM_PPP_PROTON_PHOTON_ELEMENTS);
      NewRawDiffRate(&PPPProtonElectronRate, NUC_NUM_MAIN_BINS, 
		     NUM_BG_BINS, NUM_PPP_PROTON_ELECTRON_ELEMENTS);
      NewRawDiffRate(&PPPProtonPositronRate, NUC_NUM_MAIN_BINS, 
		     NUM_BG_BINS, NUM_PPP_PROTON_POSITRON_ELEMENTS);
      NewRawDiffRate(&PPPNeutronElectronRate, NUC_NUM_MAIN_BINS, 
		     NUM_BG_BINS, NUM_PPP_PROTON_POSITRON_ELEMENTS);
      NewRawDiffRate(&PPPProtonElectronNeutrinoRate, NUC_NUM_MAIN_BINS, 
		     NUM_BG_BINS, NUM_PPP_PROTON_ANTI_MUON_NEUTRINO_ELEMENTS);
      NewRawDiffRate(&PPPProtonAntiElectronNeutrinoRate, 
		     NUC_NUM_MAIN_BINS, NUM_BG_BINS, 
		     NUM_PPP_PROTON_ANTI_ELECTRON_NEUTRINO_ELEMENTS);
      NewRawDiffRate(&PPPProtonMuonNeutrinoRate, NUC_NUM_MAIN_BINS, 
		     NUM_BG_BINS, NUM_PPP_PROTON_MUON_NEUTRINO_ELEMENTS);
      NewRawDiffRate(&PPPProtonAntiMuonNeutrinoRate, NUC_NUM_MAIN_BINS, 
		     NUM_BG_BINS, NUM_PPP_PROTON_ANTI_MUON_NEUTRINO_ELEMENTS);
      NewRawDiffRate(&PPPNeutronAntiElectronNeutrinoRate, 
		     NUC_NUM_MAIN_BINS, NUM_BG_BINS, 
		     NUM_PPP_PROTON_ANTI_MUON_NEUTRINO_ELEMENTS);
      NewRawDiffRate(&PPPNeutronMuonNeutrinoRate, NUC_NUM_MAIN_BINS, 
		     NUM_BG_BINS, NUM_PPP_PROTON_ANTI_MUON_NEUTRINO_ELEMENTS);
      NewRawDiffRate(&PPPNeutronAntiMuonNeutrinoRate, NUC_NUM_MAIN_BINS, 
		     NUM_BG_BINS, NUM_PPP_PROTON_MUON_NEUTRINO_ELEMENTS);
    }
  }

  if (NPPSwitch == 1) {
    NewRawTotalRate(&NPPTotalRate, NUC_NUM_MAIN_BINS, NUM_BG_BINS);
    
    if (nucleonToSecondarySwitch == 1)
      NewRawDiffRate(&NPPDiffRate, NUC_NUM_MAIN_BINS, NUM_BG_BINS, 
		     NUM_NPP_ELEMENTS);
  }

  // neutron decay does not fit the general category because it is not
  //    really an interaction (with background)
  if (neutronDecaySwitch == 1) {
    NewTotalRate(&neutronDecayRate, NUC_NUM_MAIN_BINS);
    NewDiffRate(&neutronDecayProtonRate, NUC_NUM_MAIN_BINS);
    
    if (nucleonToSecondarySwitch == 1)
      NewDiffRate(&neutronDecayElectronRate, NUC_NUM_MAIN_BINS);
  }

  // neutrino-neutrino rates are already folded w/ neutrino background
  if (neutrinoNeutrinoSwitch == 1) {
    NewTotalRate(&NNElNeutTotalRate, NEUT_NUM_MAIN_BINS);
    NewTotalRate(&NNMuonNeutTotalRate, NEUT_NUM_MAIN_BINS);
    NewTotalRate(&NNTauNeutTotalRate, NEUT_NUM_MAIN_BINS);
    
    NewDiffRate(&NNElNeutScatRate, NEUT_NUM_MAIN_BINS);
    NewDiffRate(&NNElNeutMuonNeutRate, NEUT_NUM_MAIN_BINS);
    NewDiffRate(&NNElNeutTauNeutRate, NEUT_NUM_MAIN_BINS);
    NewDiffRate(&NNElNeutElectronRate, NEUT_NUM_MAIN_BINS);
    NewDiffRate(&NNElNeutPhotonRate, NEUT_NUM_MAIN_BINS);
    NewDiffRate(&NNElNeutProtonRate, NEUT_NUM_MAIN_BINS);
    NewDiffRate(&NNMuonNeutScatRate, NEUT_NUM_MAIN_BINS);
    NewDiffRate(&NNMuonNeutElNeutRate, NEUT_NUM_MAIN_BINS);
    NewDiffRate(&NNMuonNeutTauNeutRate, NEUT_NUM_MAIN_BINS);
    NewDiffRate(&NNMuonNeutElectronRate, NEUT_NUM_MAIN_BINS);
    NewDiffRate(&NNMuonNeutPhotonRate, NEUT_NUM_MAIN_BINS);
    NewDiffRate(&NNMuonNeutProtonRate, NEUT_NUM_MAIN_BINS);
    NewDiffRate(&NNTauNeutScatRate, NEUT_NUM_MAIN_BINS);
    NewDiffRate(&NNTauNeutElNeutRate, NEUT_NUM_MAIN_BINS);
    NewDiffRate(&NNTauNeutMuonNeutRate, NEUT_NUM_MAIN_BINS);
    NewDiffRate(&NNTauNeutElectronRate, NEUT_NUM_MAIN_BINS);
    NewDiffRate(&NNTauNeutPhotonRate, NEUT_NUM_MAIN_BINS);
    NewDiffRate(&NNTauNeutProtonRate, NEUT_NUM_MAIN_BINS);
  }

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
  if (PPPSwitch == 1) {
    LoadPPPNucleonTables(&PPPProtonLossRate, &PPPNeutronLossRate, 
			 &PPPProtonScatRate, &PPPProtonNeutronRate, &PPPNeutronProtonRate,
			 NUM_MAIN_BINS, aDirTables);
    
    if (nucleonToSecondarySwitch == 1) {
      LoadPPPEMTables(&PPPProtonPhotonRate, &PPPProtonElectronRate, 
		      &PPPProtonPositronRate, &PPPNeutronElectronRate, 
		      NUM_MAIN_BINS, aDirTables);
      LoadPPPNeutrinoTables(&PPPProtonElectronNeutrinoRate,
			    &PPPProtonAntiElectronNeutrinoRate, &PPPProtonMuonNeutrinoRate,
			    &PPPProtonAntiMuonNeutrinoRate, 
			    &PPPNeutronAntiElectronNeutrinoRate,
			    &PPPNeutronMuonNeutrinoRate, &PPPNeutronAntiMuonNeutrinoRate, 
			    NUM_MAIN_BINS, aDirTables);
    }
  }
  if (NPPSwitch == 1) {
    LoadNPPNucleonTables(&NPPTotalRate, NUM_MAIN_BINS, aDirTables);
    
    if (nucleonToSecondarySwitch == 1)
      LoadNPPSecondaryTables(&NPPDiffRate, NUM_MAIN_BINS, aDirTables);
  }

  if (neutronDecaySwitch == 1) {
    LoadNeutronDecayNucleonTables(&neutronDecayRate, 
				  &neutronDecayProtonRate, NUM_MAIN_BINS, aDirTables);

    if (nucleonToSecondarySwitch == 1)
      LoadNeutronDecaySecondaryTables(&neutronDecayElectronRate, 
				      NUM_MAIN_BINS, aDirTables);
  }

  if (neutrinoNeutrinoSwitch == 1)
    LoadNeutrinoTables(tauNeutrinoMassSwitch, &NNElNeutTotalRate, 
		       &NNMuonNeutTotalRate, &NNTauNeutTotalRate, &NNElNeutScatRate, 
		       &NNElNeutMuonNeutRate, &NNElNeutTauNeutRate, 
		       &NNElNeutElectronRate, &NNElNeutPhotonRate, &NNElNeutProtonRate, 
		       &NNMuonNeutScatRate, &NNMuonNeutElNeutRate,
		       &NNMuonNeutTauNeutRate, &NNMuonNeutElectronRate, 
		       &NNMuonNeutPhotonRate, &NNMuonNeutProtonRate, &NNTauNeutScatRate,
		       &NNTauNeutElNeutRate, &NNTauNeutMuonNeutRate, 
		       &NNTauNeutElectronRate, &NNTauNeutPhotonRate, 
		       &NNTauNeutProtonRate, NUM_MAIN_BINS, aDirTables);
  
#ifdef DEBUG
  printf("Starting computations...\n\n");
#endif

clock_t start_comp = clock();

  //---- Initialize distance ----
  distance = 0.;  // distance variable: dist in FORTRAN
  
  PrepareSpectra(sourceTypeSwitch, &Q_0, pSpectrum, &spectrumNew, 
		 &derivative);
  // NOTE: suspect derivative is never used???

  //  if (sourceTypeSwitch == 0)  // single source
  //        DumpSpectrum(pEnergy, pEnergyWidth, &Q_0, "InitialSpectrumDump.dat"); 
  
  ComputeTotalInitialContent(pEnergy, &Q_0, &initialPhotonEnergy, 
			     &initialLeptonEnergy, &initialNucleonEnergy, 
			     &initialNeutrinoEnergy, &initialTotalEnergy, 
			     &initialPhotonNumber, &initialLeptonNumber, 
			     &initialNucleonNumber, &initialNeutrinoNumber, 
			     &initialTotalNumber);

  //--------- START of actual computation --------
  //---- initialize indices and parameters ----
  iterationCounter = 1;
  firstIndex = 0;
  leftRedshift = startingRedshift;
  deltaRedshift = 1. - (pEnergy->vector)[2]/(pEnergy->vector)[3];
  lastIndex = 0;
  propagatingDistance = 0.;
  
  
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
  
  NewDiffRate(&elNeutScatRate, NUM_MAIN_BINS);
  NewDiffRate(&elNeutMuonNeutRate, NUM_MAIN_BINS);
  NewDiffRate(&elNeutTauNeutRate, NUM_MAIN_BINS);
  NewDiffRate(&elNeutElectronRate, NUM_MAIN_BINS);
  NewDiffRate(&elNeutPhotonRate, NUM_MAIN_BINS);
  NewDiffRate(&elNeutProtonRate, NUM_MAIN_BINS);
  NewDiffRate(&muonNeutElNeutRate, NUM_MAIN_BINS);
  NewDiffRate(&muonNeutScatRate, NUM_MAIN_BINS);
  NewDiffRate(&muonNeutTauNeutRate, NUM_MAIN_BINS);
  NewDiffRate(&muonNeutElectronRate, NUM_MAIN_BINS);
  NewDiffRate(&muonNeutPhotonRate, NUM_MAIN_BINS);
  NewDiffRate(&muonNeutProtonRate, NUM_MAIN_BINS);
  NewDiffRate(&tauNeutElNeutRate, NUM_MAIN_BINS);
  NewDiffRate(&tauNeutMuonNeutRate, NUM_MAIN_BINS);
  NewDiffRate(&tauNeutScatRate, NUM_MAIN_BINS);
  NewDiffRate(&tauNeutElectronRate, NUM_MAIN_BINS);
  NewDiffRate(&tauNeutPhotonRate, NUM_MAIN_BINS);
  NewDiffRate(&tauNeutProtonRate, NUM_MAIN_BINS);
  /* I created all folded rates because I didn't want to pass all the
     switches to the subsequent functions, and these allocations are
     not so expensive anyway */
  
  //-------- This is where propagation takes place --------
  do {

    // this firstIndex & lastIndex pair is used to bin redshift finer 
    //  near the end of propagation (z close to 0) 
    ComputeRedshifts(sourceTypeSwitch, leftRedshift, &deltaRedshift,
		     &rightRedshift, &centralRedshift, &lastIndex);
    
    //---- compute various distance parameters ----
    redshiftRatio = (1. + rightRedshift)/(1. + leftRedshift);
    //    tempCoefficient = pow(OMEGA_M*pow(1. + centralRedshift,3.)+
    //			  OMEGA_LAMBDA, -1./2.)/(1.+centralRedshift);
	//        tempCoefficient = pow(1. + centralRedshift, -5./2.); 
    //    distanceStep = (leftRedshift - rightRedshift)*tempCoefficient*C/
    //      (H_0*1.e5/1.e6/DISTANCE_UNIT);

    // (cosmological parameters added July 2005)
    leftDistance = getDistance(RedshiftArray,DistanceArray,leftRedshift) ;
    rightDistance =  getDistance(RedshiftArray,DistanceArray,rightRedshift) ;
    distanceStep = leftDistance - rightDistance ;
    distance += distanceStep;
    //    rightDistance = distance - distanceStep/2.;
    //    leftDistance = distance + distanceStep/2.;
    //        rightDistance = 2./3.*C/1.e5/H_0*(1. - pow(1. + rightRedshift, -1.5));
    //        leftDistance = 2./3.*C/1.e5/H_0*(1. - pow(1. + leftRedshift, -1.5));

    evolutionFactor1 = pow((1. + centralRedshift)/(1. + startingRedshift),
			   brightPhaseExp+3.);   // redshift evolution
    numSmallSteps = (int)(ceil(distanceStep/DISTANCE_UNIT/DMAX));
    smallDistanceStep = distanceStep/DISTANCE_UNIT/numSmallSteps;
    x = 0.;
    //	printf("small stepsize: %15.6E\n", smallDistanceStep);
    
    //---- compute the photon background at given redshift ----
    LoadPhotonBackground(centralRedshift, &bgEnergy, &bgEnergyWidth, 
			 &bgPhotonDensity, aIRFlag, aZmax_IR, aRadioFlag,
			 aH0, aOmegaM, aOmegaLambda);
    
    //	if (rightRedshift < 1.e-5) 
    // DumpBgSpectrum(&bgEnergy, &bgEnergyWidth, &bgPhotonDensity,"~/");

    //---- initialize rates ----
    InitializeLeptonCoefficients(&leptonTotalRate, &leptonScatRate,
				 &leptonExchRate, &leptonPhotonRate);
    InitializePhotonCoefficients(&photonTotalRate, &photonLeptonRate);
    InitializeNucleonCoefficients(&protonTotalRate, &neutronTotalRate,
				  &protonScatRate, &protonNeutronRate, 
				  &neutronProtonRate, &protonPhotonRate, 
				  &protonElectronRate, &protonPositronRate, 
				  &neutronElectronRate, &neutronPositronRate, 
				  &protonElectronNeutrinoRate, 
				  &protonAntiElectronNeutrinoRate,
				  &protonMuonNeutrinoRate, &protonAntiMuonNeutrinoRate,
				  &neutronAntiElectronNeutrinoRate, &neutronMuonNeutrinoRate,
				  &neutronAntiMuonNeutrinoRate);
    InitializeNeutrinoCoefficients(&elNeutTotalRate, &muonNeutTotalRate,
				   &tauNeutTotalRate, &elNeutScatRate, &elNeutMuonNeutRate,
				   &elNeutTauNeutRate, &elNeutElectronRate, 
				   &elNeutPhotonRate, &elNeutProtonRate, 
				   &muonNeutElNeutRate, &muonNeutScatRate,
				   &muonNeutTauNeutRate, &muonNeutElectronRate, 
				   &muonNeutPhotonRate, &muonNeutProtonRate, 
				   &tauNeutElNeutRate, &tauNeutMuonNeutRate,
				   &tauNeutScatRate, &tauNeutElectronRate, 
				   &tauNeutPhotonRate, &tauNeutProtonRate);

    Initialize_dCVector(&continuousLoss);
    
    Initialize_dCVector(&otherLoss);
    Initialize_dCVector(&protonContinuousLoss);
    
    //---- fold interaction rates w/ photon background ----
    if (ICSSwitch == 1)
      FoldICS(&bgPhotonDensity, &ICSTotalRate, &ICSPhotonRate, 
	      &ICSScatRate, &leptonTotalRate, &leptonPhotonRate, 
	      &leptonScatRate);
    if (TPPSwitch == 1)
      FoldTPP(&bgPhotonDensity, pEnergy, &TPPTotalRate, &TPPDiffRate,
	      &leptonTotalRate, &leptonScatRate, &leptonExchRate, 
	      &otherLoss);
    if (PPSwitch == 1)
      FoldPP(&bgPhotonDensity, &PPTotalRate, &PPDiffRate, 
	     &photonTotalRate, &photonLeptonRate);
    if (DPPSwitch == 1)
      FoldDPP(&bgPhotonDensity, &DPPRate, &photonTotalRate, 
	      &photonLeptonRate);
    if (PPPSwitch == 1) {
      FoldPPPNucleon(&bgPhotonDensity, &PPPProtonLossRate,
		     &PPPNeutronLossRate, &PPPProtonScatRate, &PPPProtonNeutronRate,
		     &PPPNeutronProtonRate, &protonTotalRate, &neutronTotalRate,
		     &protonScatRate, &protonNeutronRate, &neutronProtonRate);
      
      if (nucleonToSecondarySwitch == 1)
	FoldPPPSecondary(&bgPhotonDensity, &PPPProtonPhotonRate, 
			 &PPPProtonElectronRate, &PPPProtonPositronRate, 
			 &PPPNeutronElectronRate, &PPPProtonElectronNeutrinoRate, 
			 &PPPProtonAntiElectronNeutrinoRate,
			 &PPPProtonMuonNeutrinoRate, &PPPProtonAntiMuonNeutrinoRate,
			 &PPPNeutronAntiElectronNeutrinoRate, 
			 &PPPNeutronMuonNeutrinoRate, 
			 &PPPNeutronAntiMuonNeutrinoRate, &protonPhotonRate, 
			 &protonElectronRate, &protonPositronRate, 
			 &neutronElectronRate, &neutronPositronRate, 
			 &protonElectronNeutrinoRate, 
			 &protonAntiElectronNeutrinoRate,
			 &protonMuonNeutrinoRate, &protonAntiMuonNeutrinoRate, 
			 &neutronAntiElectronNeutrinoRate, &neutronMuonNeutrinoRate,
			 &neutronAntiMuonNeutrinoRate);
    }


#ifdef EXTRACT_PHOTON_TOTAL_RATE
    // Add - E.A. Dec. 2005
    // Extract photon total energy loss rate
    for (int i=0; i<pEnergy->dimension; i++) cout << (pEnergy->vector)[i] << " " 
						  << (photonTotalRate.totalRate)[i] 
						  << endl;
    exit(0) ;
#endif
#ifdef EXTRACT_LEPTON_RATE
    for (int i=0; i<pEnergy->dimension; i++) cout << (pEnergy->vector)[i] << " " 
						  << (leptonTotalRate.totalRate)[i] 
						  << endl;
    exit(0) ;
#endif

    if (NPPSwitch == 1) {
      
      FoldNPPNucleon(&bgPhotonDensity, pEnergy, &NPPTotalRate,
		     &protonContinuousLoss);

#ifdef EXTRACT_PAIR_PROD_RATE
      // Slight add - E.A. July 2005 .
      // To extract the pair production tables of protons on photon backgrounds.
      // (need to switch NPPSwitch!)
      for (int i=0; i<pEnergy->dimension; i++) cout << (pEnergy->vector)[i] << " " 
						    << (protonContinuousLoss.vector)[i] 
						    << endl;
      exit(0) ;
#endif

      if (nucleonToSecondarySwitch == 1)
	FoldNPPSecondary(&bgPhotonDensity, &NPPDiffRate,
			 &protonPositronRate, &protonElectronRate);
    }    
	
    if (neutrinoNeutrinoSwitch == 1)
      MapNeutRates(centralRedshift, lastIndex, tauNeutrinoMassSwitch,
		   &NNElNeutTotalRate, &NNMuonNeutTotalRate, &NNTauNeutTotalRate, 
		   &NNElNeutScatRate, &NNElNeutMuonNeutRate, 
		   &NNElNeutTauNeutRate, &NNElNeutElectronRate, 
		   &NNElNeutPhotonRate, &NNElNeutProtonRate, 
		   &NNMuonNeutElNeutRate, &NNMuonNeutScatRate, 
		   &NNMuonNeutTauNeutRate, &NNMuonNeutElectronRate,
		   &NNMuonNeutPhotonRate, &NNMuonNeutProtonRate, 
		   &NNTauNeutElNeutRate, &NNTauNeutMuonNeutRate, 
		   &NNTauNeutScatRate, &NNTauNeutElectronRate,
		   &NNTauNeutPhotonRate, &NNTauNeutProtonRate, &elNeutTotalRate,
		   &muonNeutTotalRate, &tauNeutTotalRate, &elNeutScatRate,
		   &elNeutMuonNeutRate, &elNeutTauNeutRate, &elNeutElectronRate,
		   &elNeutPhotonRate, &elNeutProtonRate, &muonNeutElNeutRate,
		   &muonNeutScatRate, &muonNeutTauNeutRate, &muonNeutElectronRate,
		   &muonNeutPhotonRate, &muonNeutProtonRate, &tauNeutElNeutRate,
		   &tauNeutMuonNeutRate, &tauNeutScatRate, &tauNeutElectronRate,
		   &tauNeutPhotonRate, &tauNeutProtonRate);

    //---- main iteration (convergence) block ----
    
    for (loopCounter = 0; loopCounter < numSmallSteps; loopCounter++) {
      
      if (synchrotronSwitch == 1) {     // synchrotron == true (on)
	
	NewDiffRate(&syncRate, NUM_MAIN_BINS);
	
	B_bin = (int)((double)(B_bins)*(propagatingDistance+x)/1.e6/
		      dist_observer);
	B_loc=(pB_field->vector)[B_bin];
	InitializeSynchrotron(B_loc, pEnergy, pEnergyWidth, 
			      &synchrotronLoss, &syncRate,
			      aDirTables);
      }
      //---- compute continuous energy loss for electrons ----
      ComputeContinuousEnergyLoss(synchrotronSwitch, &synchrotronLoss, 
				  &otherLoss, &continuousLoss);
      
      evolutionFactor = evolutionFactor1;
      if ((leftDistance - x/1.e6) < minDistance)
	evolutionFactor = 0.;
      if ((leftDistance - x/1.e6) < SOURCE_CLUSTER_DISTANCE)
	evolutionFactor = SOURCE_CLUSTER_FACTOR*evolutionFactor1;
      if ((leftDistance - x/1.e6) < CLUSTER_DISTANCE)
	bkgFactor = CLUSTER_FACTOR;
      else
	bkgFactor = 1.;
      
      AdvanceNucNeutStep(sourceTypeSwitch, PPPSwitch, NPPSwitch, 
			 neutronDecaySwitch, nucleonToSecondarySwitch,
			 neutrinoNeutrinoSwitch, smallDistanceStep,
			 evolutionFactor, convergeParameter, bkgFactor, &Q_0,
			 &elNeutProtonRate,
			 &muonNeutProtonRate, &tauNeutProtonRate, &protonTotalRate,
			 &neutronTotalRate, &neutronDecayRate, &protonScatRate, 
			 &protonNeutronRate, &neutronProtonRate, 
			 &neutronDecayProtonRate, &protonMuonNeutrinoRate, 
			 &neutronAntiMuonNeutrinoRate, &protonAntiMuonNeutrinoRate, 
			 &neutronMuonNeutrinoRate, &protonElectronNeutrinoRate, 
			 &neutronAntiElectronNeutrinoRate, 
			 &protonAntiElectronNeutrinoRate, &neutronDecayElectronRate,
			 &elNeutTotalRate, &muonNeutTotalRate, &tauNeutTotalRate,
			 &elNeutScatRate, &elNeutMuonNeutRate, &elNeutTauNeutRate,
			 &muonNeutElNeutRate, &muonNeutScatRate, &muonNeutTauNeutRate,
			 &tauNeutElNeutRate, &tauNeutMuonNeutRate, &tauNeutScatRate,
			 &protonContinuousLoss, &deltaG, pSpectrum, &spectrumNew);
      
      SetNucleonSpectrum(pSpectrum, &spectrumNew);
      SetNeutrinoSpectrum(pSpectrum, &spectrumNew);
      
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
      
      if (aCutcascade_Magfield != 0 && B_loc != 0 ) {
	// Estimate the effect of B field on the 1D approximation (added E.A. June 2006)
	lEcFlag = 0 ;
	lIndex =  0;
	while (!lEcFlag) {
	  lEnergy = (pEnergy->vector)[lIndex] ;
	  // Time scales are computed in parsec
	  t_sync = 3.84e6/(lEnergy*B_loc*B_loc*ELECTRON_MASS) ;
	  t_larmor = (1.1e-21)*ELECTRON_MASS*lEnergy/(B_loc*aCutcascade_Magfield) ;
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
      iterationCounter++;
    }


    propagatingDistance += x;   // increment distance 
    
    //---- redshift bins down ----
    RedshiftDown(lastIndex, redshiftRatio, pEnergy, pSpectrum, 
		 &spectrumNew);
    
    //	printf("\nz = %g -> %g", leftRedshift, rightRedshift);
    //	printf("; d/Mpc = %g -> %g:\n", leftDistance, rightDistance);
    //	printf("%g, %g:\n", propagatingDistance, distance);
    CheckEnergy(sourceTypeSwitch, brightPhaseExp, startingRedshift,
		rightRedshift, pSpectrum, pEnergy, initialTotalEnergy);
    
    //---- prepare for new step ----
    leftRedshift = rightRedshift;
  } while (lastIndex != 1);

clock_t end_comp = clock();
  
  //---- I am done with the rates ----
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
  
  DeleteDiffRate(&elNeutScatRate);
  DeleteDiffRate(&elNeutMuonNeutRate);
  DeleteDiffRate(&elNeutTauNeutRate);
  DeleteDiffRate(&elNeutElectronRate);
  DeleteDiffRate(&elNeutPhotonRate);
  DeleteDiffRate(&elNeutProtonRate);
  DeleteDiffRate(&muonNeutElNeutRate);
  DeleteDiffRate(&muonNeutScatRate);
  DeleteDiffRate(&muonNeutTauNeutRate);
  DeleteDiffRate(&muonNeutElectronRate);
  DeleteDiffRate(&muonNeutPhotonRate);
  DeleteDiffRate(&muonNeutProtonRate);
  DeleteDiffRate(&tauNeutElNeutRate);
  DeleteDiffRate(&tauNeutMuonNeutRate);
  DeleteDiffRate(&tauNeutScatRate);
  DeleteDiffRate(&tauNeutElectronRate);
  DeleteDiffRate(&tauNeutPhotonRate);
  DeleteDiffRate(&tauNeutProtonRate);
  
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
  if (PPPSwitch == 1) {
    DeleteRawDiffRate(&PPPProtonScatRate);
    DeleteRawDiffRate(&PPPProtonNeutronRate);
    DeleteRawDiffRate(&PPPNeutronProtonRate);
    
    DeleteRawTotalRate(&PPPProtonLossRate);
    DeleteRawTotalRate(&PPPNeutronLossRate);
    
    if (nucleonToSecondarySwitch == 1) {
      DeleteRawDiffRate(&PPPProtonPhotonRate);
      DeleteRawDiffRate(&PPPProtonElectronRate);
      DeleteRawDiffRate(&PPPProtonPositronRate);
      DeleteRawDiffRate(&PPPNeutronElectronRate);
      DeleteRawDiffRate(&PPPProtonElectronNeutrinoRate);
      DeleteRawDiffRate(&PPPProtonAntiElectronNeutrinoRate);
      DeleteRawDiffRate(&PPPProtonMuonNeutrinoRate);
      DeleteRawDiffRate(&PPPProtonAntiMuonNeutrinoRate);
      DeleteRawDiffRate(&PPPNeutronAntiElectronNeutrinoRate);
      DeleteRawDiffRate(&PPPNeutronMuonNeutrinoRate);
      DeleteRawDiffRate(&PPPNeutronAntiMuonNeutrinoRate);
    }
  }
  if (NPPSwitch == 1) {
    if (nucleonToSecondarySwitch == 1)
      DeleteRawDiffRate(&NPPDiffRate);
    DeleteRawTotalRate(&NPPTotalRate);
  }
  if (DPPSwitch == 1)
    DeleteRawTotalRate(&DPPRate);
  
  if (neutronDecaySwitch == 1) {
    if (nucleonToSecondarySwitch == 1)
      DeleteDiffRate(&neutronDecayElectronRate);
    DeleteDiffRate(&neutronDecayProtonRate);
    DeleteTotalRate(&neutronDecayRate);
  }
  
  
  if (neutrinoNeutrinoSwitch == 1) {
    DeleteDiffRate(&NNElNeutScatRate);
    DeleteDiffRate(&NNElNeutMuonNeutRate);
    DeleteDiffRate(&NNElNeutTauNeutRate);
    DeleteDiffRate(&NNElNeutElectronRate);
    DeleteDiffRate(&NNElNeutPhotonRate);
    DeleteDiffRate(&NNElNeutProtonRate);
    DeleteDiffRate(&NNMuonNeutElNeutRate);
    DeleteDiffRate(&NNMuonNeutScatRate);
    DeleteDiffRate(&NNMuonNeutTauNeutRate);
    DeleteDiffRate(&NNMuonNeutElectronRate);
    DeleteDiffRate(&NNMuonNeutPhotonRate);
    DeleteDiffRate(&NNMuonNeutProtonRate);
    DeleteDiffRate(&NNTauNeutElNeutRate);
    DeleteDiffRate(&NNTauNeutMuonNeutRate);
    DeleteDiffRate(&NNTauNeutScatRate);
    DeleteDiffRate(&NNTauNeutElectronRate);
    DeleteDiffRate(&NNTauNeutPhotonRate);
    DeleteDiffRate(&NNTauNeutProtonRate);
    
    DeleteTotalRate(&NNElNeutTotalRate);
    DeleteTotalRate(&NNMuonNeutTotalRate);
    DeleteTotalRate(&NNTauNeutTotalRate);
  }
  
  //    FinalPrintOutToTheScreen(distance, startingRedshift, propagatingDistance);
  
  DeleteSpectrum(&Q_0);
  DeleteSpectrum(&spectrumNew);
  DeleteSpectrum(&derivative);
  
  Delete_dCVector(&synchrotronLoss);
  Delete_dCVector(&continuousLoss);
  Delete_dCVector(&otherLoss);
  Delete_dCVector(&protonContinuousLoss);
  
  Delete_dCVector(&deltaG);
  Delete_dCVector(&bgEnergy);
  Delete_dCVector(&bgEnergyWidth);
  Delete_dCVector(&bgPhotonDensity);

  Delete_dCVector(&RedshiftArray) ;
  Delete_dCVector(&DistanceArray) ;
clock_t end = clock();

clock_t total = end -start;
clock_t comp = end_comp -start_comp;

#ifdef DEBUG
printf("Dint Comptime: %f\n", float(comp)/float(total));
#endif
}

void BuildRedshiftTable(double aH0, double aOmegaM, double aOmegaLambda, 
			dCVector* pRedshiftArray, dCVector* pDistanceArray) {

  // Routine added Jan 2006 - E.A.

  // aH0 in km/s/Mpc
  // C is in cm/s, distances in Mpc
  int lNbBins = NBINS_REDSHIFTTABLE ;
  double lZmin = ZMIN_REDSHIFTTABLE ;
  double lZmax = ZMAX_REDSHIFTTABLE ;
  New_dCVector(pRedshiftArray, lNbBins) ;
  New_dCVector(pDistanceArray, lNbBins) ;

  pRedshiftArray->vector[0] = 0 ;
  for (int i=0; i < lNbBins - 1; i++) {
    double kk = i/(lNbBins-2.) ;
    pRedshiftArray->vector[i+1] = lZmin*pow(lZmax/lZmin,kk) ;
  }
  dCVector lpEvolutionFactor ;
  New_dCVector(&lpEvolutionFactor, lNbBins) ;
  lpEvolutionFactor.vector[0] = 1 ;
  pDistanceArray->vector[0] = 0 ;
  for (int i=1; i < lNbBins; i++) {
    double lZ = pRedshiftArray->vector[i] ;
    lpEvolutionFactor.vector[i] = 1./((1.+lZ)*sqrt(aOmegaLambda+aOmegaM*pow(1.+lZ,3))) ;
    pDistanceArray->vector[i] = pDistanceArray->vector[i-1]+
      (lZ-pRedshiftArray->vector[i-1])*0.5*
      (lpEvolutionFactor.vector[i-1]+lpEvolutionFactor.vector[i])*C/(aH0*1.e5) ;
  }
  Delete_dCVector(&lpEvolutionFactor) ;

}

double getRedshift(dCVector RedshiftArray, dCVector DistanceArray, double distance) {

  // Routine added Jan 2006 - E.A.

  unsigned int i0 = 0 ;
  while (DistanceArray.vector[i0+1] <= distance) i0 += 1 ;
  double lRedshift = RedshiftArray.vector[i0] +  (RedshiftArray.vector[i0+1]-RedshiftArray.vector[i0])
    *(distance-DistanceArray.vector[i0])
    /(DistanceArray.vector[i0+1]-DistanceArray.vector[i0]) ;
  if ( distance <= 0 ) lRedshift = 0 ;

  return lRedshift ;
}

double getDistance(dCVector RedshiftArray, dCVector DistanceArray, double redshift) {

  // Routine added Jan 2006 - E.A.

  unsigned int i0 = 0 ;
  while (RedshiftArray.vector[i0+1] <= redshift) i0 += 1 ;
  double lDistance = DistanceArray.vector[i0] +  (DistanceArray.vector[i0+1]-DistanceArray.vector[i0])
    *(redshift-RedshiftArray.vector[i0])
    /(RedshiftArray.vector[i0+1]-RedshiftArray.vector[i0]) ;
  if ( redshift <= 0 ) lDistance = 0 ;

  // At that level lDistance is in Mpc --> we convert it into cm.
  lDistance *= 1.e6 ; // in pc
  lDistance *= DISTANCE_UNIT ; // in cm

  return lDistance ;
}



#include "dint/advance.h"

#ifdef DEBUG
void PrintEnergy(const dCVector* pArray, const dCVector* pEnergy)
{
    double temp = 0.;
    int i;

    if (pArray->dimension != pEnergy->dimension)
    {
	Error("PrintEnergy: inconsistent dimensions", PROGRAM_ERROR);
    }

    for (i = 0; i < pEnergy->dimension; i++)
    {
	temp += (pArray->vector)[i]*(pEnergy->vector)[i];
    }

    printf("%15.6E\n", temp);
}
#endif

void ComputeRedshifts(const int sourceTypeSwitch, const double leftRedshift,
		      double* pDeltaRedshift, double* pRightRedshift,
		      double* pCentralRedshift, int* pLastIndex) {
    double clusterRedshift,localRedshift,minRedshift,maxRedshift;

    clusterRedshift = pow(1.-1.5e5*H_0/C_0*CLUSTER_DISTANCE,-2./3.)-1.;
    localRedshift = pow(1.-1.5e5*H_0/C_0*SOURCE_CLUSTER_DISTANCE,-2./3.)-1.;
    maxRedshift = localRedshift;
    minRedshift = maxRedshift;
    if (clusterRedshift > maxRedshift)
      maxRedshift = clusterRedshift;
    else
      minRedshift = clusterRedshift;
    if (sourceTypeSwitch == 0)
      minRedshift = maxRedshift;

    if ((*pLastIndex == 0) && (leftRedshift < 0.6)) {
	*pDeltaRedshift /= 2.;
	*pLastIndex = 3;
    }
    *pRightRedshift = leftRedshift - (1. + leftRedshift)*(*pDeltaRedshift);
    if (*pRightRedshift < 0.) {
	*pRightRedshift = 0.;
	*pLastIndex = 1;
    }
    *pCentralRedshift = (leftRedshift + *pRightRedshift)/2.;
    if ((*pLastIndex == 1) && (*pCentralRedshift > maxRedshift)) {
	*pRightRedshift = maxRedshift;
	*pCentralRedshift = (leftRedshift + *pRightRedshift)/2.;
	*pLastIndex = 2;
    }

    if ((*pLastIndex == 1) && (*pCentralRedshift > minRedshift)) {
	*pRightRedshift = minRedshift;
	*pCentralRedshift = (leftRedshift + *pRightRedshift)/2.;
	*pLastIndex = 2;
    }

}

void AdvanceNucleonStep(const int sourceTypeSwitch,
			const int PPPSwitch, const int NPPSwitch,
			const int neutronDecaySwitch,
			const int neutrinoNeutrinoSwitch,
			const double smallDistanceStep,
			const double evolutionFactor,
			const double convergeParameter,
			const double bkgFactor,
			const Spectrum* pQ_0,
			const DiffRate* elNeutProtonRate,
			const DiffRate* muonNeutProtonRate,
			const DiffRate* tauNeutProtonRate,
			const TotalRate* protonTotalRate,
			const TotalRate* neutronTotalRate,
			const TotalRate* neutronDecayRate,
			const DiffRate* protonScatRate,
			const DiffRate* protonNeutronRate,
			const DiffRate* neutronProtonRate,
			const DiffRate* neutronDecayProtonRate,
			const dCVector* protonContinuousLoss,
			const dCVector* deltaG, const Spectrum* pSpectrum,
			Spectrum* pSpectrumNew) {
  Spectrum spectrumTemp;
  Spectrum influx;
  Spectrum influx0;
  Spectrum influxExt;
  Spectrum outflux;
  int num_main_bins;
  double changeMax;

  num_main_bins = pSpectrum->numberOfMainBins;

  NewSpectrum(&spectrumTemp, num_main_bins);
  NewSpectrum(&influx, num_main_bins);
  NewSpectrum(&influx0, num_main_bins);
  NewSpectrum(&influxExt, num_main_bins);
  NewSpectrum(&outflux, num_main_bins);

  GetExternalFlux(sourceTypeSwitch, evolutionFactor, PROTON, pQ_0,
		  &influxExt);
  GetExternalFlux(sourceTypeSwitch, evolutionFactor, NEUTRON, pQ_0,
		  &influxExt);

  if (neutrinoNeutrinoSwitch == 1)
    GetNucleonInfluxFromNeutrinos(bkgFactor, elNeutProtonRate,
				  muonNeutProtonRate, tauNeutProtonRate,
				  pSpectrumNew, &influx0);

  do {
    SetNucleonSpectrum(&spectrumTemp, pSpectrumNew);

    InitializeSpectrum(&influx);
    InitializeSpectrum(&outflux);

    if ((PPPSwitch == 1) || (NPPSwitch == 1) ||
	(neutronDecaySwitch == 1)) {
      // call all the processes that produce nucleons from nucleons
      GetNucleonFluxFromNucleons(neutronDecaySwitch, protonTotalRate,
				 neutronTotalRate, neutronDecayRate, protonScatRate,
				 neutronProtonRate, protonNeutronRate,
				 neutronDecayProtonRate, protonContinuousLoss,
				 deltaG, pSpectrumNew, &influx, &outflux);
    }

    ImplicitEquation(smallDistanceStep, PROTON, &influx, &influx0,
		     &influxExt, &outflux, pSpectrum, pSpectrumNew);
    ImplicitEquation(smallDistanceStep, NEUTRON, &influx, &influx0,
		     &influxExt, &outflux, pSpectrum, pSpectrumNew);

    changeMax = 0.;
    ComputeChange(&spectrumTemp, pSpectrumNew, PROTON, &changeMax);
    ComputeChange(&spectrumTemp, pSpectrumNew, NEUTRON, &changeMax);
  } while (changeMax > convergeParameter);

  DeleteSpectrum(&spectrumTemp);
  DeleteSpectrum(&influx);
  DeleteSpectrum(&influx0);
  DeleteSpectrum(&influxExt);
  DeleteSpectrum(&outflux);
}


void AdvanceNeutrinoStep(const int sourceTypeSwitch,
			 const int neutrinoNeutrinoSwitch,
			 const int PPPSwitch,
			 const int neutronDecaySwitch,
			 const int nucleonToSecondarySwitch,
			 const double smallDistanceStep,
			 const double evolutionFactor,
			 const double convergeParameter,
			 const double bkgFactor,
			 const Spectrum* pQ_0,
			 const DiffRate* protonMuonNeutrinoRate,
			 const DiffRate* neutronAntiMuonNeutrinoRate,
			 const DiffRate* protonAntiMuonNeutrinoRate,
			 const DiffRate* neutronMuonNeutrinoRate,
			 const DiffRate* protonElectronNeutrinoRate,
			 const DiffRate* neutronAntiElectronNeutrinoRate,
			 const DiffRate* protonAntiElectronNeutrinoRate,
			 const DiffRate* neutronDecayElectronRate,
			 const TotalRate* elNeutTotalRate,
			 const TotalRate* muonNeutTotalRate,
			 const TotalRate* tauNeutTotalRate,
			 const DiffRate* elNeutScatRate,
			 const DiffRate* elNeutMuonNeutRate,
			 const DiffRate* elNeutTauNeutRate,
			 const DiffRate* muonNeutElNeutRate,
			 const DiffRate* muonNeutScatRate,
			 const DiffRate* muonNeutTauNeutRate,
			 const DiffRate* tauNeutElNeutRate,
			 const DiffRate* tauNeutMuonNeutRate,
			 const DiffRate* tauNeutScatRate,
			 const Spectrum* pSpectrum, Spectrum* pSpectrumNew) {
    double changeMax;
    int num_main_bins;
    Spectrum spectrumTemp;
    Spectrum influx;
    Spectrum influx0;
    Spectrum influxExt;
    Spectrum outflux;

    num_main_bins = pSpectrum->numberOfMainBins;

    NewSpectrum(&spectrumTemp, num_main_bins);
    NewSpectrum(&influx, num_main_bins);
    NewSpectrum(&influx0, num_main_bins);
    NewSpectrum(&influxExt, num_main_bins);
    NewSpectrum(&outflux, num_main_bins);

    //    SetNeutrinoSpectrum(spectrumNew, spectrum);

    GetExternalFlux(sourceTypeSwitch, evolutionFactor, MUON_NEUTRINO, pQ_0,
        &influxExt);
    GetExternalFlux(sourceTypeSwitch, evolutionFactor, ANTI_MUON_NEUTRINO,
        pQ_0, &influxExt);
    GetExternalFlux(sourceTypeSwitch, evolutionFactor, ELECTRON_NEUTRINO, pQ_0,
        &influxExt);
    GetExternalFlux(sourceTypeSwitch, evolutionFactor, ANTI_ELECTRON_NEUTRINO,
        pQ_0, &influxExt);
    GetExternalFlux(sourceTypeSwitch, evolutionFactor, TAU_NEUTRINO, pQ_0,
        &influxExt);
    GetExternalFlux(sourceTypeSwitch, evolutionFactor, ANTI_TAU_NEUTRINO, pQ_0,
        &influxExt);

    InitializeSpectrum(&influx0);
    /* this is not really necessary but will keep it for redundancy */

    /* neutrino influx from nucleons */
    if (nucleonToSecondarySwitch == 1)
    {
        if ((PPPSwitch == 1) || (neutronDecaySwitch == 1))
	{
	    GetNeutrinoInfluxFromNucleons(neutronDecaySwitch,
                protonMuonNeutrinoRate, neutronAntiMuonNeutrinoRate,
                neutronMuonNeutrinoRate, protonAntiMuonNeutrinoRate,
                protonElectronNeutrinoRate, neutronAntiElectronNeutrinoRate,
                protonAntiElectronNeutrinoRate, neutronDecayElectronRate,
                pSpectrumNew, &influx0);
	}
    }

    do {
	SetNeutrinoSpectrum(&spectrumTemp, pSpectrumNew);

	InitializeSpectrum(&influx);
	InitializeSpectrum(&outflux);

	/* call processes that produce neutrinos from neutrinos */
	if (neutrinoNeutrinoSwitch == 1)
	{
	    GetNeutrinoFluxFromNeutrinos(bkgFactor, elNeutTotalRate,
		muonNeutTotalRate,
                tauNeutTotalRate, elNeutScatRate, elNeutMuonNeutRate,
                elNeutTauNeutRate, muonNeutElNeutRate, muonNeutScatRate,
                muonNeutTauNeutRate, tauNeutElNeutRate, tauNeutMuonNeutRate,
                tauNeutScatRate, pSpectrumNew, &influx, &outflux);
	}

	ImplicitEquation(smallDistanceStep, MUON_NEUTRINO, &influx,
            &influx0, &influxExt, &outflux, pSpectrum, pSpectrumNew);
	ImplicitEquation(smallDistanceStep, ANTI_MUON_NEUTRINO, &influx,
            &influx0, &influxExt, &outflux, pSpectrum, pSpectrumNew);
	ImplicitEquation(smallDistanceStep, ELECTRON_NEUTRINO, &influx,
            &influx0, &influxExt, &outflux, pSpectrum, pSpectrumNew);
	ImplicitEquation(smallDistanceStep, ANTI_ELECTRON_NEUTRINO, &influx,
            &influx0, &influxExt, &outflux, pSpectrum, pSpectrumNew);
	ImplicitEquation(smallDistanceStep, TAU_NEUTRINO, &influx,
            &influx0, &influxExt, &outflux, pSpectrum, pSpectrumNew);
	ImplicitEquation(smallDistanceStep, ANTI_TAU_NEUTRINO, &influx,
            &influx0, &influxExt, &outflux, pSpectrum, pSpectrumNew);

	changeMax = 0.;
	ComputeChange(&spectrumTemp, pSpectrumNew, ELECTRON_NEUTRINO,
            &changeMax);
	ComputeChange(&spectrumTemp, pSpectrumNew, ANTI_ELECTRON_NEUTRINO,
            &changeMax);
	ComputeChange(&spectrumTemp, pSpectrumNew, MUON_NEUTRINO,
            &changeMax);
	ComputeChange(&spectrumTemp, pSpectrumNew, ANTI_MUON_NEUTRINO,
            &changeMax);
	ComputeChange(&spectrumTemp, pSpectrumNew, TAU_NEUTRINO,
            &changeMax);
	ComputeChange(&spectrumTemp, pSpectrumNew, ANTI_TAU_NEUTRINO,
            &changeMax);
    } while (changeMax > convergeParameter);

    DeleteSpectrum(&spectrumTemp);
    DeleteSpectrum(&influx);
    DeleteSpectrum(&influx0);
    DeleteSpectrum(&influxExt);
    DeleteSpectrum(&outflux);
}

void AdvanceNucNeutStep(const int sourceTypeSwitch,
			const int PPPSwitch, const int NPPSwitch,
			const int neutronDecaySwitch,
			const int nucleonToSecondarySwitch,
			const int neutrinoNeutrinoSwitch,
			const double smallDistanceStep,
			const double evolutionFactor,
			const double convergeParameter,
			const double bkgFactor,
			const Spectrum* pQ_0,
			const DiffRate* elNeutProtonRate,
			const DiffRate* muonNeutProtonRate,
			const DiffRate* tauNeutProtonRate,
			const TotalRate* protonTotalRate,
			const TotalRate* neutronTotalRate,
			const TotalRate* neutronDecayRate,
			const DiffRate* protonScatRate,
			const DiffRate* protonNeutronRate,
			const DiffRate* neutronProtonRate,
			const DiffRate* neutronDecayProtonRate,
			const DiffRate* protonMuonNeutrinoRate,
			const DiffRate* neutronAntiMuonNeutrinoRate,
			const DiffRate* protonAntiMuonNeutrinoRate,
			const DiffRate* neutronMuonNeutrinoRate,
			const DiffRate* protonElectronNeutrinoRate,
			const DiffRate* neutronAntiElectronNeutrinoRate,
			const DiffRate* protonAntiElectronNeutrinoRate,
			const DiffRate* neutronDecayElectronRate,
			const TotalRate* elNeutTotalRate,
			const TotalRate* muonNeutTotalRate,
			const TotalRate* tauNeutTotalRate,
			const DiffRate* elNeutScatRate,
			const DiffRate* elNeutMuonNeutRate,
			const DiffRate* elNeutTauNeutRate,
			const DiffRate* muonNeutElNeutRate,
			const DiffRate* muonNeutScatRate,
			const DiffRate* muonNeutTauNeutRate,
			const DiffRate* tauNeutElNeutRate,
			const DiffRate* tauNeutMuonNeutRate,
			const DiffRate* tauNeutScatRate,
			const dCVector* protonContinuousLoss,
			const dCVector* deltaG, const Spectrum* pSpectrum,
			Spectrum* pSpectrumNew)
{
    Spectrum spectrumGlobalTemp;
    /*  temporary spectrum used in swapping spectra */
    double changeGlobal;
    int num_main_bins;

    num_main_bins = pSpectrum->numberOfMainBins;

    NewSpectrum(&spectrumGlobalTemp, num_main_bins);

    do {
	SetNucleonSpectrum(&spectrumGlobalTemp, pSpectrumNew);
	SetNeutrinoSpectrum(&spectrumGlobalTemp, pSpectrumNew);


	AdvanceNucleonStep(sourceTypeSwitch, PPPSwitch, NPPSwitch,
            neutronDecaySwitch, neutrinoNeutrinoSwitch, smallDistanceStep,
            evolutionFactor, convergeParameter, bkgFactor, pQ_0,
            elNeutProtonRate, muonNeutProtonRate, tauNeutProtonRate,
            protonTotalRate, neutronTotalRate, neutronDecayRate,
            protonScatRate, protonNeutronRate, neutronProtonRate,
            neutronDecayProtonRate, protonContinuousLoss, deltaG, pSpectrum,
            pSpectrumNew);

	AdvanceNeutrinoStep(sourceTypeSwitch, neutrinoNeutrinoSwitch,
            PPPSwitch, neutronDecaySwitch, nucleonToSecondarySwitch,
            smallDistanceStep, evolutionFactor, convergeParameter, bkgFactor,
            pQ_0, protonMuonNeutrinoRate, neutronAntiMuonNeutrinoRate,
            protonAntiMuonNeutrinoRate, neutronMuonNeutrinoRate,
            protonElectronNeutrinoRate, neutronAntiElectronNeutrinoRate,
            protonAntiElectronNeutrinoRate, neutronDecayElectronRate,
            elNeutTotalRate, muonNeutTotalRate, tauNeutTotalRate,
            elNeutScatRate, elNeutMuonNeutRate, elNeutTauNeutRate,
            muonNeutElNeutRate, muonNeutScatRate, muonNeutTauNeutRate,
            tauNeutElNeutRate, tauNeutMuonNeutRate, tauNeutScatRate,
            pSpectrum, pSpectrumNew);

        /*global convergence */
        changeGlobal = 0.;
	ComputeChange(&spectrumGlobalTemp, pSpectrumNew, PROTON,
            &changeGlobal);
	ComputeChange(&spectrumGlobalTemp, pSpectrumNew, NEUTRON,
            &changeGlobal);
	ComputeChange(&spectrumGlobalTemp, pSpectrumNew, ELECTRON_NEUTRINO,
            &changeGlobal);
	ComputeChange(&spectrumGlobalTemp, pSpectrumNew,
            ANTI_ELECTRON_NEUTRINO, &changeGlobal);
	ComputeChange(&spectrumGlobalTemp, pSpectrumNew, MUON_NEUTRINO,
            &changeGlobal);
	ComputeChange(&spectrumGlobalTemp, pSpectrumNew, ANTI_MUON_NEUTRINO,
            &changeGlobal);
	ComputeChange(&spectrumGlobalTemp, pSpectrumNew, TAU_NEUTRINO,
            &changeGlobal);
	ComputeChange(&spectrumGlobalTemp, pSpectrumNew, ANTI_TAU_NEUTRINO,
            &changeGlobal);
    } while (changeGlobal > convergeParameter);


    DeleteSpectrum(&spectrumGlobalTemp);
}


void AdvanceEMStep(const int sourceTypeSwitch, const int PPSwitch,
		   const int ICSSwitch, const int TPPSwitch,
		   const int DPPSwitch, const int synchrotronSwitch,
		   const int PPPSwitch, const int NPPSwitch,
		   const int neutronDecaySwitch,
		   const int nucleonToSecondarySwitch,
		   const int neutrinoNeutrinoSwitch,
		   const double smallDistanceStep,
		   const double evolutionFactor,
		   const double convergeParameter, const double bkgFactor,
		   const Spectrum* pQ_0,
		   const DiffRate* photonLeptonRate,
		   const DiffRate* protonElectronRate,
		   const DiffRate* neutronPositronRate,
		   const DiffRate* protonPositronRate,
		   const DiffRate* neutronElectronRate,
		   const DiffRate* neutronDecayElectronRate,
		   const DiffRate* elNeutElectronRate,
		   const DiffRate* muonNeutElectronRate,
		   const DiffRate* tauNeutElectronRate,
		   const DiffRate* protonPhotonRate,
		   const DiffRate* elNeutPhotonRate,
		   const DiffRate* muonNeutPhotonRate,
		   const DiffRate* tauNeutPhotonRate,
		   const TotalRate* leptonTotalRate,
		   const DiffRate* leptonScatRate,
		   const DiffRate* leptonExchRate,
		   const dCVector* continuousLoss, const dCVector* deltaG,
		   const TotalRate* photonTotalRate,
		   const DiffRate* leptonPhotonRate,
		   const DiffRate* syncRate, const Spectrum* pSpectrum,
		   Spectrum* pSpectrumNew)
{
    int i;
    Spectrum spectrumTemp;
    Spectrum spectrumGlobalTemp;
    /* temporary spectra used in swapping spectra */
    Spectrum influx;
    Spectrum influx0;
    Spectrum influxExt;
    Spectrum outflux;
    double changeMax;
    /* maximum difference by advance: diff in FORTRAN */
    double changeGlobal;
    int num_main_bins;

    num_main_bins = pSpectrum->numberOfMainBins;

    NewSpectrum(&spectrumTemp, num_main_bins);
    NewSpectrum(&spectrumGlobalTemp, num_main_bins);
    NewSpectrum(&influx, num_main_bins);
    NewSpectrum(&influx0, num_main_bins);
    NewSpectrum(&influxExt, num_main_bins);
    NewSpectrum(&outflux, num_main_bins);

    do {
        SetEMSpectrum(&spectrumGlobalTemp, pSpectrumNew);

	//        for (i = 0; i < NUM_MAIN_BINS; i++)
	//            spectrumOld[PHOTON][i] = spectrumNew[PHOTON][i];
	/* I don't think I need it here... */

        /*---- Converge e+/- first ----*/
	InitializeSpectrum(&influxExt);

        /* external injection */
	GetExternalFlux(sourceTypeSwitch, evolutionFactor, ELECTRON, pQ_0,
	    &influxExt);
	GetExternalFlux(sourceTypeSwitch, evolutionFactor, POSITRON, pQ_0,
	    &influxExt);


	InitializeSpectrum(&influx0);

        /* influx from pair production/double pair production */
	if ((PPSwitch == 1) || (DPPSwitch == 1))
	{
	    GetLeptonInfluxFromPhotons(photonLeptonRate, pSpectrumNew,
                &influx0);
	}

	/* influx from nucleons */
	if (nucleonToSecondarySwitch == 1)    /* secondary tables included */
	{
	    if ((PPPSwitch == 1) || (NPPSwitch == 1) ||
		(neutronDecaySwitch == 1))
	    {
	        GetLeptonInfluxFromNucleons(neutronDecaySwitch,
                    protonElectronRate, neutronPositronRate,
                    protonPositronRate, neutronElectronRate,
                    neutronDecayElectronRate, pSpectrumNew, &influx0);
	    }
	}

	/* influx from neutrinos */
	if (neutrinoNeutrinoSwitch == 1)
	{
	    GetLeptonInfluxFromNeutrinos(bkgFactor, elNeutElectronRate,
                muonNeutElectronRate, tauNeutElectronRate, pSpectrumNew,
                &influx0);
	}

        do {
            for (i = 0; i < num_main_bins; i++)
            {
                spectrumTemp.spectrum[ELECTRON][i] =
		    (pSpectrumNew->spectrum)[ELECTRON][i];
                spectrumTemp.spectrum[POSITRON][i] =
		    (pSpectrumNew->spectrum)[POSITRON][i];
            }

	    InitializeSpectrum(&influx);
	    InitializeSpectrum(&outflux);

            /* influx & outflux from inverse Compton scattering/synchrotron
                radiation */
	    if ((ICSSwitch == 1) || (TPPSwitch == 1) ||
		(synchrotronSwitch == 1))
	    {
		GetLeptonFluxFromLeptons(leptonTotalRate, leptonScatRate,
                    leptonExchRate, continuousLoss, deltaG, pSpectrumNew,
		    &influx, &outflux);
	    }

	    ImplicitEquation(smallDistanceStep, ELECTRON, &influx, &influx0,
		&influxExt, &outflux, pSpectrum, pSpectrumNew);
	    ImplicitEquation(smallDistanceStep, POSITRON, &influx, &influx0,
		&influxExt, &outflux, pSpectrum, pSpectrumNew);
            /* main difference equation part */

            changeMax = 0.;
	    ComputeChange(&spectrumTemp, pSpectrumNew, ELECTRON, &changeMax);
	    ComputeChange(&spectrumTemp, pSpectrumNew, POSITRON, &changeMax);
        } while (changeMax > convergeParameter);


        /*---- Converge photons ----*/
        /* external influx */
	InitializeSpectrum(&influxExt);

	GetExternalFlux(sourceTypeSwitch, evolutionFactor, PHOTON, pQ_0,
	    &influxExt);

	InitializeSpectrum(&influx0);

        /* influx from electrons (ICS/synchrotron) */
	if ((ICSSwitch == 1) || (synchrotronSwitch == 1))
	{
	    GetPhotonInfluxFromLeptons(leptonPhotonRate, synchrotronSwitch,
                syncRate, pSpectrumNew, &influx0);
	}

	/* influx from nucleons */
	if ((nucleonToSecondarySwitch == 1) && (PPPSwitch == 1))
	{
	    GetPhotonInfluxFromNucleons(protonPhotonRate, pSpectrumNew,
                &influx0);
	}
	if (neutrinoNeutrinoSwitch == 1)
	{
	    GetPhotonInfluxFromNeutrinos(bkgFactor, elNeutPhotonRate,
		muonNeutPhotonRate, tauNeutPhotonRate, pSpectrumNew, &influx0);
	}

        do {
            for (i = 0; i < num_main_bins; i++)
	    {
                spectrumTemp.spectrum[PHOTON][i] =
		    (pSpectrumNew->spectrum)[PHOTON][i];
	    }

            /* loss (and zero gain) by PP/DPP */
	    InitializeSpectrum(&outflux);
	    if ((PPSwitch == 1) || (DPPSwitch == 1))
	    {
		GetPhotonFluxFromPhotons(photonTotalRate, &outflux);
	    }

	    ImplicitEquation(smallDistanceStep, PHOTON, &influx, &influx0,
		&influxExt, &outflux, pSpectrum, pSpectrumNew);
            /* main difference equation */

            changeMax = 0.;     /* reset it for photons */
	    ComputeChange(&spectrumTemp, pSpectrumNew, PHOTON, &changeMax);
        } while (changeMax > convergeParameter);

        /*global convergence */
        changeGlobal = 0.;
	ComputeChange(&spectrumGlobalTemp, pSpectrumNew, ELECTRON,
            &changeGlobal);
	ComputeChange(&spectrumGlobalTemp, pSpectrumNew, POSITRON,
            &changeGlobal);
	ComputeChange(&spectrumGlobalTemp, pSpectrumNew, PHOTON,
            &changeGlobal);
    } while (changeGlobal > convergeParameter);

    DeleteSpectrum(&spectrumTemp);
    DeleteSpectrum(&spectrumGlobalTemp);
    DeleteSpectrum(&influx);
    DeleteSpectrum(&influx0);
    DeleteSpectrum(&influxExt);
    DeleteSpectrum(&outflux);
}


void RedshiftDown(const int lastIndex, const double redshiftRatio,
                  const dCVector* pEnergy, Spectrum* pSpectrum,
		  Spectrum* pSpectrumNew)
{
    int i;

    if (pEnergy->dimension != pSpectrum->numberOfMainBins)
    {
	Error("RedshiftDown: inconsistent dimensions", PROGRAM_ERROR);
    }
    for (i = 0; i < NUM_SPECIES; i++)
    {
	RedshiftBinsDown(lastIndex, redshiftRatio, pEnergy,
	    pSpectrum->spectrum[i], pSpectrumNew->spectrum[i]);
    }
}


void RedshiftBinsDown(const int lastIndex, const double redshiftRatio,
                      const dCVector* pEnergy, double* pSpectrum,
		      double* pSpectrumNew)
/* redshift particle distributions at the end of a big dx */
{
    int i;
    int num_main_bins;
    double averaging_factor;

    averaging_factor = (pow(10., 1./2./BINS_PER_DECADE) +
        pow(10., -1./2./BINS_PER_DECADE))/2.;

    num_main_bins = pEnergy->dimension;

    if (lastIndex == 0) /* normal steps far from end */
    {
        for (i = 0; i < num_main_bins-1; i++)
        {
            pSpectrum[i] = pSpectrum[i+1]*redshiftRatio*
                redshiftRatio*redshiftRatio;
        }
        pSpectrum[num_main_bins-1] = 0.;
    }
    else
    /* if it's close to z = 0, need some fine tuning because redshift
        steps are taken more finely */
    {
        //    double spectrumTemp[NUM_MAIN_BINS];
        dCVector spectrumTemp;
        double energyTemp;
        double inflow;
        int iMod;
        int iMod2;
        double tempFraction;
        double tempExponent;

	New_dCVector(&spectrumTemp, num_main_bins);

        for (i = 0; i < num_main_bins; i++)
        {
            energyTemp = (pEnergy->vector)[i]/redshiftRatio;
            if (energyTemp > (pEnergy->vector)[num_main_bins-1]/
		averaging_factor*pow(10., 1./2./BINS_PER_DECADE))
	    {
                inflow = 0.;
	    }
            else
            {
		iMod = (int)(BINS_PER_DECADE*(log10(energyTemp*ELECTRON_MASS) -
		    MAX_ENERGY_EXP) + (num_main_bins + 0.5));
#ifdef DEBUG
		CheckIndex(0, num_main_bins, iMod, "RedshiftBinsDown");
#endif
                if ((pEnergy->vector)[iMod] > energyTemp)
		{
                    iMod--;
		}
                if (iMod < 0)
		{
                    Error("RedshiftBinsDown: iMod became negative.",
			PROGRAM_ERROR);   /* get the hell out of here! */
		}
                iMod2 = iMod + 1;
                if (iMod2 >= num_main_bins)
		{
                    inflow = pSpectrum[num_main_bins-1];
		}
                else if ((pSpectrum[iMod] < 1.e-40) ||
                    (pSpectrum[iMod2] < 1.e-40))
                {
#ifdef DEBUG
		    CheckIndex(0, num_main_bins, iMod2, "RedshiftBinsDown");
#endif
                    tempFraction = (energyTemp - (pEnergy->vector)[iMod])/
                        ((pEnergy->vector)[iMod2] - (pEnergy->vector)[iMod]);
                    if (tempFraction < 0.)
		    {
                        tempFraction = 0.;
		    }
                    inflow = pSpectrum[iMod]*(1. - tempFraction) +
                        pSpectrum[iMod2]*tempFraction;
                }
                else
                {
                    tempExponent = log(energyTemp/(pEnergy->vector)[iMod])/
                        log((pEnergy->vector)[iMod2]/(pEnergy->vector)[iMod]);
                    inflow = pSpectrum[iMod]*pow(pSpectrum[iMod2]/
                        pSpectrum[iMod], tempExponent);
                }
            }
            spectrumTemp.vector[i] = inflow*redshiftRatio*redshiftRatio*
                redshiftRatio;
        }

#ifdef DEBUG
	if (spectrumTemp.vector[num_main_bins-1] > pSpectrum[num_main_bins-1])
	    printf("Auch!!!\n");
#endif

        for (i = 0; i < num_main_bins; i++)
        {
	    if (spectrumTemp.vector[i] < 0.)
	    {
		Error("RedshiftBinsDown: spectrumTemp negative.",
		    PROGRAM_ERROR);
	    }
            pSpectrum[i] = spectrumTemp.vector[i];
            pSpectrumNew[i] = pSpectrum[i];
	}
	Delete_dCVector(&spectrumTemp);
    }
}


void GetExternalFlux(const int sourceTypeSwitch, const double evolutionFactor,
		     const PARTICLE particle, const Spectrum* pQ_0,
		     Spectrum* pInfluxExt)
{
    int i;

    if (pQ_0->numberOfMainBins != pInfluxExt->numberOfMainBins)
    {
	Error("GetExternalFlux: inconsistent dimensions", PROGRAM_ERROR);
    }

    for (i = 0; i < pQ_0->numberOfMainBins; i++)
    {
	if (sourceTypeSwitch == 0)	/* single source */
	{
	    (pInfluxExt->spectrum)[particle][i] = 0.;
	}
	else
	{
	    (pInfluxExt->spectrum)[particle][i] =
	        (pQ_0->spectrum)[particle][i]*evolutionFactor;
	}
    }
}

void ImplicitEquation(const double smallDistanceStep,
		      const PARTICLE particle, const Spectrum* pInflux,
		      const Spectrum* pInflux0, const Spectrum* pInfluxExt,
		      const Spectrum* pOutflux, const Spectrum* pSpectrum,
		      Spectrum* pSpectrumNew)
{
    /* influx: influx from the same KIND of species (e-: e+/-, p: p/n, etc.)
       influx0: influx from other species
       influxExt: external influx
       NOTE: there is no fundamental distinction between influx and influx0 */
    int i;
    int num_main_bins;

    num_main_bins = pSpectrum->numberOfMainBins;

    if ((pInflux->numberOfMainBins != num_main_bins) ||
	(pInflux0->numberOfMainBins != num_main_bins) ||
	(pInfluxExt->numberOfMainBins != num_main_bins) ||
	(pOutflux->numberOfMainBins != num_main_bins) ||
	(pSpectrumNew->numberOfMainBins != num_main_bins))
    {
	Error("ImplicitEquation: inconsistent dimensions", PROGRAM_ERROR);
    }

    for (i = 0; i < num_main_bins; i++)
    {
	(pSpectrumNew->spectrum)[particle][i] =
	    ((pInflux->spectrum)[particle][i] +
            (pInflux0->spectrum)[particle][i] +
	    (pInfluxExt->spectrum)[particle][i] +
            (pSpectrum->spectrum)[particle][i]/smallDistanceStep)/
	    ((pOutflux->spectrum)[particle][i] + 1./smallDistanceStep);
	/* main difference equation */

	if ((pSpectrumNew->spectrum)[particle][i] < 0.)
	{
	    (pSpectrumNew->spectrum)[particle][i] = 0.;
	}
	/* reset if it passes below zero */
    }
}

void ExplicitEquation(const double smallDistanceStep,
		      const PARTICLE particle, const Spectrum* pInflux,
		      const Spectrum* pInflux0, const Spectrum* pInfluxExt,
		      const Spectrum* pOutflux, const Spectrum* pSpectrum,
		      const Spectrum* pSpectrumNew)
{
    int i;
    int num_main_bins;

    num_main_bins = pSpectrum->numberOfMainBins;

    if ((pInflux->numberOfMainBins != num_main_bins) ||
	(pInflux0->numberOfMainBins != num_main_bins) ||
	(pInfluxExt->numberOfMainBins != num_main_bins) ||
	(pOutflux->numberOfMainBins != num_main_bins) ||
	(pSpectrumNew->numberOfMainBins != num_main_bins))
    {
	Error("ImplicitEquation: inconsistent dimensions", PROGRAM_ERROR);
    }

    for (i = 0; i < num_main_bins; i++)
    {
	(pSpectrumNew->spectrum)[particle][i] =
	    (pSpectrum->spectrum)[particle][i] + smallDistanceStep*
	    ((pInflux->spectrum)[particle][i] +
	    (pInflux0->spectrum)[particle][i] +
            (pInfluxExt->spectrum)[particle][i] -
	    (pSpectrum->spectrum)[particle][i]*
            (pOutflux->spectrum)[particle][i]);

	if ((pSpectrumNew->spectrum)[particle][i] < 0.)
	{
	    (pSpectrumNew->spectrum)[particle][i] = 0.;
	}
    }
}

void ComputeChange(const Spectrum* pSpectrumTemp,
		   const Spectrum* pSpectrumNew,
		   const PARTICLE particle, double* pChangeMax)
{
    double change;
    int i;

    if (pSpectrumTemp->numberOfMainBins != pSpectrumNew->numberOfMainBins)
    {
	Error("ComputeChange: inconsistent dimensions", PROGRAM_ERROR);
    }

    for (i = 0; i < pSpectrumNew->numberOfMainBins; i++)
    {
	change = fabs((pSpectrumNew->spectrum)[particle][i] -
	    (pSpectrumTemp->spectrum)[particle][i])/
	    (fabs((pSpectrumTemp->spectrum)[particle][i]) + EPSILON);
	if (change > *pChangeMax)
	{
	    *pChangeMax = change;
	}
    }
}

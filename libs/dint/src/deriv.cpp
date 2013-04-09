#include "dint/deriv.h"
#include "dint/rate.h"
#include "dint/spectrum.h"
#include "dint/error.h"
#include "dint/check.h"
#include "dint/cvector.h"


void GetLeptonInfluxFromPhotons(const DiffRate* photonLeptonRate,
				const Spectrum* pSpectrumNew, 
				Spectrum* pInflux0)
{
    ComputeInflux(1., photonLeptonRate, PHOTON, ELECTRON, pSpectrumNew, 
		  pInflux0);
    ComputeInflux(1., photonLeptonRate, PHOTON, POSITRON, pSpectrumNew, 
		  pInflux0);
}

void GetLeptonInfluxFromNucleons(const int neutronDecaySwitch,
				 const DiffRate* protonElectronRate,
				 const DiffRate* neutronPositronRate,
				 const DiffRate* protonPositronRate,
				 const DiffRate* neutronElectronRate,
				 const DiffRate* neutronDecayElectronRate,
				 const Spectrum* pSpectrumNew,
				 Spectrum* pInflux0)
{
    /* lepton influx from nucleons */
    ComputeInflux(1., protonElectronRate, PROTON, ELECTRON, pSpectrumNew, 
	pInflux0);
    ComputeInflux(1., neutronPositronRate, NEUTRON, POSITRON, pSpectrumNew, 
	pInflux0);
    ComputeInflux(1., protonPositronRate, PROTON, POSITRON, pSpectrumNew, 
	pInflux0);
    ComputeInflux(1., neutronElectronRate, NEUTRON, ELECTRON, pSpectrumNew, 
	pInflux0);
    if (neutronDecaySwitch == 1)
    {
	ComputeInflux(1., neutronDecayElectronRate, NEUTRON, ELECTRON, 
            pSpectrumNew, pInflux0);
    }
}

void GetLeptonInfluxFromNeutrinos(const double bkgFactor,
				  const DiffRate* elNeutElectronRate,
				  const DiffRate* muonNeutElectronRate,
				  const DiffRate* tauNeutElectronRate,
				  const Spectrum* pSpectrumNew,
				  Spectrum* pInflux0)
{
    ComputeInflux(bkgFactor, elNeutElectronRate, ELECTRON_NEUTRINO, ELECTRON, 
        pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, elNeutElectronRate, ANTI_ELECTRON_NEUTRINO,
	ELECTRON, pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, elNeutElectronRate, ELECTRON_NEUTRINO, POSITRON, 
        pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, elNeutElectronRate, ANTI_ELECTRON_NEUTRINO,
	POSITRON, pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, muonNeutElectronRate, MUON_NEUTRINO, ELECTRON, 
        pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, muonNeutElectronRate, ANTI_MUON_NEUTRINO,
        ELECTRON, pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, muonNeutElectronRate, MUON_NEUTRINO, POSITRON, 
        pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, muonNeutElectronRate, ANTI_MUON_NEUTRINO,
        POSITRON, pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, tauNeutElectronRate, TAU_NEUTRINO, ELECTRON, 
        pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, tauNeutElectronRate, ANTI_TAU_NEUTRINO, ELECTRON, 
        pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, tauNeutElectronRate, TAU_NEUTRINO, POSITRON, 
        pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, tauNeutElectronRate, ANTI_TAU_NEUTRINO, POSITRON, 
        pSpectrumNew, pInflux0);
}

void GetLeptonFluxFromLeptons(const TotalRate* leptonTotalRate,
			      const DiffRate* leptonScatRate,
			      const DiffRate* leptonExchRate,
			      const dCVector* continuousLoss, 
			      const dCVector* deltaG,
			      const Spectrum* pSpectrumNew, Spectrum* pInflux, 
			      Spectrum* pOutflux)
{
    int i;
    int num_main_bins;

    num_main_bins = pOutflux->numberOfMainBins;
    if ((continuousLoss->dimension != num_main_bins) ||
	(deltaG->dimension != num_main_bins))
    {
	Error("GetLeptonFluxFromLeptons: inconsistent dimensions", 
	    PROGRAM_ERROR);
    }

    /* ICS & TPP (differential) */
    ComputeOutflux(1., leptonTotalRate, ELECTRON, pOutflux);
    ComputeOutflux(1., leptonTotalRate, POSITRON, pOutflux);
    /* outflux */

    ComputeInflux(1., leptonScatRate, ELECTRON, ELECTRON, pSpectrumNew,
        pInflux);
    ComputeInflux(1., leptonScatRate, POSITRON, POSITRON, pSpectrumNew,
        pInflux);
    ComputeInflux(1., leptonExchRate, ELECTRON, POSITRON, pSpectrumNew,
        pInflux);
    ComputeInflux(1., leptonExchRate, POSITRON, ELECTRON, pSpectrumNew,
        pInflux);
    /* this is for TPP */

    /* continuous energy loss (synchrotron & TPP) */
    for (i = 0; i < num_main_bins; i++)
    {
        double fraction;

        fraction = (continuousLoss->vector)[i]/(deltaG->vector)[i];
        (pOutflux->spectrum)[ELECTRON][i] += -fraction;
        (pOutflux->spectrum)[POSITRON][i] += -fraction;
        if (i != num_main_bins - 1)
        {
            (pInflux->spectrum)[ELECTRON][i] += 
	        -(continuousLoss->vector)[i+1]/(deltaG->vector)[i+1]*
		(pSpectrumNew->spectrum)[ELECTRON][i+1];
            (pInflux->spectrum)[POSITRON][i] += 
	        -(continuousLoss->vector)[i+1]/(deltaG->vector)[i+1]*
		(pSpectrumNew->spectrum)[POSITRON][i+1];
        }
    }
}


void GetPhotonInfluxFromLeptons(const DiffRate* leptonPhotonRate,
				const int synchrotronSwitch,
				const DiffRate* syncRate,
				const Spectrum* pSpectrumNew, 
				Spectrum* pInflux0)
{
    int i;
    int num_main_bins;

    num_main_bins = pInflux0->numberOfMainBins;
    if (synchrotronSwitch == 1)
    { 
	if (syncRate->mainDimension != num_main_bins)
	{
	    Error("GetPhotonInfluxFromLeptons: inconsistent dimensions",
	          PROGRAM_ERROR);
	}
    }

    ComputeInflux(1., leptonPhotonRate, ELECTRON, PHOTON, pSpectrumNew,
        pInflux0);
    ComputeInflux(1., leptonPhotonRate, POSITRON, PHOTON, pSpectrumNew,
        pInflux0);
    /* influx from ICS */

    /* influx from synchrotron */
    if (synchrotronSwitch == 1)
    {
        for (i = 0; i < num_main_bins; i++)
        {
	    if (((syncRate->bound)[i][0] != num_main_bins - 1) ||
		(syncRate->bound)[i][1] != 0)
            {
                int j;

                for (j = (syncRate->bound)[i][0]; j <= (syncRate->bound)[i][1];
		     j++)
                {
                    (pInflux0->spectrum)[PHOTON][j] += 
		        (syncRate->diffRate)[i][j]*
                        ((pSpectrumNew->spectrum)[ELECTRON][i] + 
			(pSpectrumNew->spectrum)[POSITRON][i]);
                }
            }
        }
    }
}

void GetPhotonInfluxFromNucleons(const DiffRate* protonPhotonRate,
				 const Spectrum* pSpectrumNew, 
				 Spectrum* pInflux0)
{
    ComputeInflux(1., protonPhotonRate, PROTON, PHOTON, pSpectrumNew, 
        pInflux0);
    ComputeInflux(1., protonPhotonRate, NEUTRON, PHOTON, pSpectrumNew, 
        pInflux0);
}

void GetPhotonInfluxFromNeutrinos(const double bkgFactor,
				  const DiffRate* elNeutPhotonRate,
				  const DiffRate* muonNeutPhotonRate,
				  const DiffRate* tauNeutPhotonRate,
				  const Spectrum* pSpectrumNew,
				  Spectrum* pInflux0)
{
    ComputeInflux(bkgFactor, elNeutPhotonRate, ELECTRON_NEUTRINO, PHOTON,
	pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, elNeutPhotonRate, ANTI_ELECTRON_NEUTRINO, PHOTON, 
        pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, muonNeutPhotonRate, MUON_NEUTRINO, PHOTON,
        pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, muonNeutPhotonRate, ANTI_MUON_NEUTRINO, PHOTON, 
        pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, tauNeutPhotonRate, TAU_NEUTRINO, PHOTON,
        pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, tauNeutPhotonRate, ANTI_TAU_NEUTRINO, PHOTON, 
        pSpectrumNew, pInflux0);
}


void GetPhotonFluxFromPhotons(const TotalRate* photonTotalRate, 
			      Spectrum* pOutflux)
{
    ComputeOutflux(1., photonTotalRate, PHOTON, pOutflux);
}

void GetNucleonFluxFromNucleons(const int neutronDecaySwitch,
				const TotalRate* protonTotalRate, 
				const TotalRate* neutronTotalRate,
				const TotalRate* neutronDecayRate,
				const DiffRate* protonScatRate,
				const DiffRate* neutronProtonRate,
				const DiffRate* protonNeutronRate,
				const DiffRate* neutronDecayProtonRate,
				const dCVector* protonContinuousLoss,
				const dCVector* deltaG, 
				const Spectrum* pSpectrumNew,
				Spectrum* pInflux, Spectrum* pOutflux)
{
    int i;
    int num_main_bins;

    num_main_bins = pInflux->numberOfMainBins;
    if ((protonContinuousLoss->dimension != num_main_bins) ||
	(deltaG->dimension != num_main_bins))
    {
	Error("GetNucleonFluxFromNucleons: inconsistent dimensions",
	      PROGRAM_ERROR);
    }

    ComputeOutflux(1., protonTotalRate, PROTON, pOutflux);
    ComputeOutflux(1., neutronTotalRate, NEUTRON, pOutflux);
    if (neutronDecaySwitch == 1)
    {
	ComputeOutflux(1., neutronDecayRate, NEUTRON, pOutflux);
    }
    /* compute outflux */

    ComputeInflux(1., protonScatRate, PROTON, PROTON, pSpectrumNew, pInflux);
    ComputeInflux(1., neutronProtonRate, NEUTRON, PROTON, pSpectrumNew,
        pInflux);
    ComputeInflux(1., protonScatRate, NEUTRON, NEUTRON, pSpectrumNew, pInflux);
    ComputeInflux(1., protonNeutronRate, PROTON, NEUTRON, pSpectrumNew,
        pInflux);
    /* compute influx by PPP */

    if (neutronDecaySwitch == 1)
    {
      /*ComputeInflux(1., neutronDecayProtonRate, NEUTRON, PROTON, pSpectrumNew,
	pInflux);*/

    for (i = 0; i < neutronDecayRate->mainDimension; i++)
    {
	pInflux->spectrum[PROTON][i] += (neutronDecayRate->totalRate[i])*
        (pSpectrumNew->spectrum)[NEUTRON][i];
    }

    }
    /* add proton influx by neutron decay */

    for (i = 0; i < num_main_bins; i++)
    {
	double fraction;

	fraction = (protonContinuousLoss->vector)[i]/(deltaG->vector)[i];
	(pOutflux->spectrum)[PROTON][i] += -fraction;
	if (i != num_main_bins - 1)
	{
	    (pInflux->spectrum)[PROTON][i] += -fraction*
	        (pSpectrumNew->spectrum)[PROTON][i+1];
	}
    }
}

void GetNeutrinoInfluxFromNucleons(const int neutronDecaySwitch,
				   const DiffRate* protonMuonNeutrinoRate,
				   const DiffRate* neutronAntiMuonNeutrinoRate,
				   const DiffRate* neutronMuonNeutrinoRate,
				   const DiffRate* protonAntiMuonNeutrinoRate,
				   const DiffRate* protonElectronNeutrinoRate,
				   const DiffRate* neutronAntiElectronNeutrinoRate,
				   const DiffRate* protonAntiElectronNeutrinoRate,
				   const DiffRate* neutronDecayElectronRate,
				   const Spectrum* pSpectrumNew,
				   Spectrum* pInflux0)
{
    ComputeInflux(1., protonMuonNeutrinoRate, PROTON, MUON_NEUTRINO,
        pSpectrumNew, pInflux0);
    ComputeInflux(1., neutronAntiMuonNeutrinoRate, NEUTRON,
        ANTI_MUON_NEUTRINO, pSpectrumNew, pInflux0);
    ComputeInflux(1., protonAntiMuonNeutrinoRate, PROTON, ANTI_MUON_NEUTRINO, 
        pSpectrumNew, pInflux0);
    ComputeInflux(1., neutronMuonNeutrinoRate, NEUTRON, MUON_NEUTRINO, 
        pSpectrumNew, pInflux0);
    ComputeInflux(1., protonElectronNeutrinoRate, PROTON, ELECTRON_NEUTRINO, 
        pSpectrumNew, pInflux0);
    ComputeInflux(1., neutronAntiElectronNeutrinoRate, NEUTRON, 
        ANTI_ELECTRON_NEUTRINO, pSpectrumNew, pInflux0);
    ComputeInflux(1., protonAntiElectronNeutrinoRate, PROTON, 
        ANTI_ELECTRON_NEUTRINO, pSpectrumNew, pInflux0);
    ComputeInflux(1., protonAntiElectronNeutrinoRate, NEUTRON,
        ELECTRON_NEUTRINO, pSpectrumNew, pInflux0);

    if (neutronDecaySwitch == 1)
    {
	ComputeInflux(1., neutronDecayElectronRate, NEUTRON, 
            ANTI_ELECTRON_NEUTRINO, pSpectrumNew, pInflux0);
    }
}

void GetNucleonInfluxFromNeutrinos(const double bkgFactor,
				   const DiffRate* elNeutProtonRate,
				   const DiffRate* muonNeutProtonRate,
				   const DiffRate* tauNeutProtonRate,
				   const Spectrum* pSpectrumNew, 
				   Spectrum* pInflux0)
{
    ComputeInflux(bkgFactor, elNeutProtonRate, ELECTRON_NEUTRINO, PROTON,
        pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, elNeutProtonRate, ANTI_ELECTRON_NEUTRINO, PROTON, 
        pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, muonNeutProtonRate, MUON_NEUTRINO, PROTON,
        pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, muonNeutProtonRate, ANTI_MUON_NEUTRINO, PROTON,
        pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, tauNeutProtonRate, TAU_NEUTRINO, PROTON,
        pSpectrumNew, pInflux0);
    ComputeInflux(bkgFactor, tauNeutProtonRate, ANTI_TAU_NEUTRINO, PROTON,
        pSpectrumNew, pInflux0);

    /* neutron is unchanged: this has to be replaced after proton-antiproton 
       symmetry is solved */
}

void GetNeutrinoFluxFromNeutrinos(const double bkgFactor,
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
				  const Spectrum* pSpectrumNew,
				  Spectrum* pInflux, Spectrum* pOutflux)
{
    int i;
    int k;
    int num_main_bins;

    num_main_bins = pInflux->numberOfMainBins;

    ComputeOutflux(bkgFactor, elNeutTotalRate, ELECTRON_NEUTRINO, pOutflux);
    ComputeOutflux(bkgFactor, elNeutTotalRate, ANTI_ELECTRON_NEUTRINO,
        pOutflux);
    ComputeOutflux(bkgFactor, muonNeutTotalRate, MUON_NEUTRINO, pOutflux);
    ComputeOutflux(bkgFactor, muonNeutTotalRate, ANTI_MUON_NEUTRINO, pOutflux);
    ComputeOutflux(bkgFactor, tauNeutTotalRate, TAU_NEUTRINO, pOutflux);
    ComputeOutflux(bkgFactor, tauNeutTotalRate, ANTI_TAU_NEUTRINO, pOutflux);


    ComputeInflux(bkgFactor, elNeutScatRate, ELECTRON_NEUTRINO,
        ELECTRON_NEUTRINO, pSpectrumNew, pInflux);
    ComputeInflux(bkgFactor, elNeutScatRate, ANTI_ELECTRON_NEUTRINO,
        ELECTRON_NEUTRINO, pSpectrumNew, pInflux);
    ComputeInflux(bkgFactor, muonNeutElNeutRate, MUON_NEUTRINO,
        ELECTRON_NEUTRINO, pSpectrumNew, pInflux);
    ComputeInflux(bkgFactor, muonNeutElNeutRate, ANTI_MUON_NEUTRINO,
        ELECTRON_NEUTRINO, pSpectrumNew, pInflux);
    ComputeInflux(bkgFactor, tauNeutElNeutRate, TAU_NEUTRINO,
        ELECTRON_NEUTRINO, pSpectrumNew, pInflux);
    ComputeInflux(bkgFactor, tauNeutElNeutRate, ANTI_TAU_NEUTRINO,
        ELECTRON_NEUTRINO, pSpectrumNew, pInflux);

    ComputeInflux(bkgFactor, elNeutMuonNeutRate, ELECTRON_NEUTRINO,
        MUON_NEUTRINO, pSpectrumNew, pInflux);
    ComputeInflux(bkgFactor, elNeutMuonNeutRate, ANTI_ELECTRON_NEUTRINO,
        MUON_NEUTRINO, pSpectrumNew, pInflux);
    ComputeInflux(bkgFactor, muonNeutScatRate, MUON_NEUTRINO, MUON_NEUTRINO,
        pSpectrumNew, pInflux);
    ComputeInflux(bkgFactor, muonNeutScatRate, ANTI_MUON_NEUTRINO,
        MUON_NEUTRINO, pSpectrumNew, pInflux);
    ComputeInflux(bkgFactor, tauNeutMuonNeutRate, TAU_NEUTRINO, MUON_NEUTRINO,
        pSpectrumNew, pInflux);
    ComputeInflux(bkgFactor, tauNeutMuonNeutRate, ANTI_TAU_NEUTRINO,
        MUON_NEUTRINO, pSpectrumNew, pInflux);

    ComputeInflux(bkgFactor, elNeutTauNeutRate, ELECTRON_NEUTRINO,
        TAU_NEUTRINO, pSpectrumNew, pInflux);
    ComputeInflux(bkgFactor, elNeutTauNeutRate, ANTI_ELECTRON_NEUTRINO,
        TAU_NEUTRINO, pSpectrumNew, pInflux);
    ComputeInflux(bkgFactor, muonNeutTauNeutRate, MUON_NEUTRINO, TAU_NEUTRINO,
        pSpectrumNew, pInflux);
    ComputeInflux(bkgFactor, muonNeutTauNeutRate, ANTI_MUON_NEUTRINO,
        TAU_NEUTRINO, pSpectrumNew, pInflux);
    ComputeInflux(bkgFactor, tauNeutScatRate, TAU_NEUTRINO, TAU_NEUTRINO,
        pSpectrumNew, pInflux);
    ComputeInflux(bkgFactor, tauNeutScatRate, ANTI_TAU_NEUTRINO, TAU_NEUTRINO,
        pSpectrumNew, pInflux);

    /* ??? */
    for (i = 0; i < num_main_bins; i++)
    {
	for (k = 0; k <= i; k++)
	{
	    (pInflux->spectrum)[ANTI_ELECTRON_NEUTRINO][k] = 
	        (pInflux->spectrum)[ELECTRON_NEUTRINO][k];
	    (pInflux->spectrum)[ANTI_MUON_NEUTRINO][k] = 
	        (pInflux->spectrum)[MUON_NEUTRINO][k];
	    (pInflux->spectrum)[ANTI_TAU_NEUTRINO][k] = 
	        (pInflux->spectrum)[TAU_NEUTRINO][k];
	}
    }
}

void ComputeOutflux(const double bkgFactor, const TotalRate* pRate,
		    const PARTICLE parent, Spectrum* pOutflux)
{
    int i;

    if (pRate->mainDimension != pOutflux->numberOfMainBins)
    {
	Error("ComputeOutflux: inconsistent dimension", PROGRAM_ERROR);
    }

    for (i = 0; i < pRate->mainDimension; i++)
    {
	pOutflux->spectrum[parent][i] += bkgFactor*(pRate->totalRate[i]);
    }
}

void ComputeInflux(const double bkgFactor, const DiffRate* pRate,
		   const PARTICLE parent, const PARTICLE daughter,
		   const Spectrum* pSpectrum, Spectrum* pInflux)
{
    int i;
    int k;

    if ((pRate->mainDimension != pSpectrum->numberOfMainBins) ||
	(pSpectrum->numberOfMainBins != pInflux->numberOfMainBins))
    {
	Error("ComputeInflux: inconsistent dimensions", PROGRAM_ERROR);
    }

    for (i = 0; i < pRate->mainDimension; i++)
    {
	if ((pRate->bound[i][0] != pRate->mainDimension - 1) ||
	    (pRate->bound[i][1] != 0))
	/* bound is valid */
	{
	    for (k = pRate->bound[i][0]; k <= pRate->bound[i][1]; k++)
	    {
#ifdef DEBUG
		CheckIndex(0, i+1, k, "ComputeInflux");
#endif
		pInflux->spectrum[daughter][k] += bkgFactor*
		   (pRate->diffRate)[i][k]*
		   (pSpectrum->spectrum)[parent][i];
	    }
	}
    }
}

#include <stdarg.h>
#include <math.h>
#include "dint/rate.h"
#include "dint/const.h"
#include "dint/cvector.h"
#include "dint/check.h"
#include "dint/utilities.h"


void InitializeLeptonCoefficients(TotalRate* leptonTotalRate,
				  DiffRate* leptonScatRate,
                                  DiffRate* leptonExchRate,
                                  DiffRate* leptonPhotonRate)
{
    InitializeTotalRate(leptonTotalRate);
    InitializeDiffRate(leptonScatRate);
    InitializeDiffRate(leptonExchRate);
    InitializeDiffRate(leptonPhotonRate);
}


void InitializePhotonCoefficients(TotalRate* photonTotalRate,
                                  DiffRate* photonLeptonRate)
{
    InitializeTotalRate(photonTotalRate);
    InitializeDiffRate(photonLeptonRate);
}

void InitializeNucleonCoefficients(TotalRate* protonTotalRate,
				   TotalRate* neutronTotalRate,
				   DiffRate* protonScatRate,
				   DiffRate* protonNeutronRate,
				   DiffRate* neutronProtonRate,
				   DiffRate* protonPhotonRate,
				   DiffRate* protonElectronRate,
				   DiffRate* protonPositronRate,
				   DiffRate* neutronElectronRate,
				   DiffRate* neutronPositronRate,
				   DiffRate* protonElectronNeutrinoRate,
				   DiffRate* protonAntiElectronNeutrinoRate,
				   DiffRate* protonMuonNeutrinoRate,
				   DiffRate* protonAntiMuonNeutrinoRate,
				   DiffRate* neutronAntiElectronNeutrinoRate,
				   DiffRate* neutronMuonNeutrinoRate,
				   DiffRate* neutronAntiMuonNeutrinoRate)
{
    InitializeTotalRate(protonTotalRate);
    InitializeTotalRate(neutronTotalRate);

    InitializeDiffRate(protonScatRate);
    InitializeDiffRate(protonNeutronRate);
    InitializeDiffRate(neutronProtonRate);
    InitializeDiffRate(protonPhotonRate);
    InitializeDiffRate(protonElectronRate);
    InitializeDiffRate(protonPositronRate);
    InitializeDiffRate(neutronElectronRate);
    InitializeDiffRate(neutronPositronRate);
    InitializeDiffRate(protonElectronNeutrinoRate);
    InitializeDiffRate(protonAntiElectronNeutrinoRate);
    InitializeDiffRate(protonMuonNeutrinoRate);
    InitializeDiffRate(protonAntiMuonNeutrinoRate);
    InitializeDiffRate(neutronAntiElectronNeutrinoRate);
    InitializeDiffRate(neutronMuonNeutrinoRate);
    InitializeDiffRate(neutronAntiMuonNeutrinoRate);
}

void InitializeNeutrinoCoefficients(TotalRate* elNeutTotalRate,
				    TotalRate* muonNeutTotalRate,
				    TotalRate* tauNeutTotalRate,
				    DiffRate* elNeutScatRate,
				    DiffRate* elNeutMuonNeutRate,
				    DiffRate* elNeutTauNeutRate,
				    DiffRate* elNeutElectronRate,
				    DiffRate* elNeutPhotonRate,
				    DiffRate* elNeutProtonRate,
				    DiffRate* muonNeutElNeutRate,
				    DiffRate* muonNeutScatRate,
				    DiffRate* muonNeutTauNeutRate,
				    DiffRate* muonNeutElectronRate,
				    DiffRate* muonNeutPhotonRate,
				    DiffRate* muonNeutProtonRate,
				    DiffRate* tauNeutElNeutRate,
				    DiffRate* tauNeutMuonNeutRate,
				    DiffRate* tauNeutScatRate,
				    DiffRate* tauNeutElectronRate,
				    DiffRate* tauNeutPhotonRate,
				    DiffRate* tauNeutProtonRate)
{
    InitializeTotalRate(elNeutTotalRate);
    InitializeTotalRate(muonNeutTotalRate);
    InitializeTotalRate(tauNeutTotalRate);

    InitializeDiffRate(elNeutScatRate);
    InitializeDiffRate(elNeutMuonNeutRate);
    InitializeDiffRate(elNeutTauNeutRate);
    InitializeDiffRate(elNeutElectronRate);
    InitializeDiffRate(elNeutPhotonRate);
    InitializeDiffRate(elNeutProtonRate);
    InitializeDiffRate(muonNeutElNeutRate);
    InitializeDiffRate(muonNeutScatRate);
    InitializeDiffRate(muonNeutTauNeutRate);
    InitializeDiffRate(muonNeutElectronRate);
    InitializeDiffRate(muonNeutPhotonRate);
    InitializeDiffRate(muonNeutProtonRate);
    InitializeDiffRate(tauNeutElNeutRate);
    InitializeDiffRate(tauNeutMuonNeutRate);
    InitializeDiffRate(tauNeutScatRate);
    InitializeDiffRate(tauNeutElectronRate);
    InitializeDiffRate(tauNeutPhotonRate);
    InitializeDiffRate(tauNeutProtonRate);
}


void FoldTotalRate(const dCVector* pBgPhotonDensity, 
		   const RawTotalRate* pRawTotalRate, TotalRate* pTotalRate)
{
    int i;
    int j;

    if ((pBgPhotonDensity->dimension != pRawTotalRate->bgDimension) ||
	(pRawTotalRate->mainDimension != pTotalRate->mainDimension))
    {
	Error("FoldTotalRate: inconsistent dimensions", PROGRAM_ERROR);
    }

    for (i = 0; i < pRawTotalRate->mainDimension; i++)
    {
	for (j = 0; j < pRawTotalRate->bgDimension; j++)
	{
	    (pTotalRate->totalRate)[i] += (pBgPhotonDensity->vector)[j]*
	        (pRawTotalRate->totalRate)[i][j];
	}
    }
}

void FoldDiffRate(const dCVector* pBgPhotonDensity,
		  const RawDiffRate* pRawDiffRate,
		  DiffRate* pDiffRate, const int scatSwitch, ...)
{
    int i;
    int j;
    int k;
    int offset = 0;
    int jLower;
    int jUpper;
    int num_main_bins;
    int num_bg_bins;

    num_main_bins = pRawDiffRate->mainDimension;
    num_bg_bins = pRawDiffRate->bgDimension;

    if ((pBgPhotonDensity->dimension != num_bg_bins) ||
	(num_main_bins != pDiffRate->mainDimension))
    {
        Error("FoldDiffRate: inconsistent dimensions", PROGRAM_ERROR);
    }

    if (scatSwitch != 0)
    /* scattering type: adjustment of total rate(s) needed */
    {
	va_list pArg;
//	TotalRate* totalRateArray[scatSwitch];
	TotalRate* totalRateArray[3];

	va_start(pArg, scatSwitch);
	for (i = 0; i < scatSwitch; i++)
	    totalRateArray[i] = va_arg(pArg, TotalRate*);
	/* Let totalRateArray[i] point to the totalRate* we will modify; we do
	   not need to update the arguments later
	   (i.e. va_arg(...) = totalRateArray[i];) because the actual values
	   they point to have been properly updated */

	for (i = 0; i < num_main_bins; i++)
	{
	    jLower = num_main_bins - 1;
	    jUpper = 0;
	    for (j = 0; j < num_bg_bins; j++)
	    {
		if (pRawDiffRate->bound[i][j][0] != -1)
		/* above threshold */
		{
		    offset += -(pRawDiffRate->bound)[i][j][0];
		    jLower = IMin(jLower, (pRawDiffRate->bound)[i][j][0]);
		    jUpper = IMax(jUpper, (pRawDiffRate->bound)[i][j][1]);
		    for (k = (pRawDiffRate->bound)[i][j][0]; 
			 k <= (pRawDiffRate->bound)[i][j][1]; k++)
		    {
#ifdef DEBUG
			CheckIndex(0, i+1, k, "FoldDiffRate");
			CheckIndex(0, pRawDiffRate->numberOfElements, k+offset,
				   "FoldDiffRate");
#endif
			(pDiffRate->diffRate)[i][k] += 
			    (pBgPhotonDensity->vector)[j]*
			    (pRawDiffRate->diffRate)[k+offset];
		    }
		    /* these lines take care of the appropriate subtractions 
		       made in the main implicit formula */
		    if ((pRawDiffRate->bound)[i][j][1] == i)
		    /*  if (((pRawDiffRate->bound)[i][j][0] <= i) && 
			    ((pRawDiffRate->bound)[i][j][1] >= i))    */
		    {
			int l;

			for (l = 0; l < scatSwitch; l++)
			{
			    (totalRateArray[l]->totalRate)[i] += 
			        -(pRawDiffRate->diffRate)[i+offset]*
			        (pBgPhotonDensity->vector)[j];
			}
			(pDiffRate->diffRate)[i][i] += 
			    -(pRawDiffRate->diffRate)[i+offset]*
			    (pBgPhotonDensity->vector)[j];
		    }
		    offset += (pRawDiffRate->bound)[i][j][1] + 1;
		    /* reset the offset index */
		}
	    }
	    (pDiffRate->bound)[i][0] = IMin((pDiffRate->bound)[i][0], jLower);
	    (pDiffRate->bound)[i][1] = IMax((pDiffRate->bound)[i][1], jUpper);
        }
	va_end(pArg);
    }
    else
    {
	for (i = 0; i < num_main_bins; i++)
	{
	    jLower = num_main_bins - 1;
	    jUpper = 0;
	    for (j = 0; j < num_bg_bins; j++)
	    {
		if (pRawDiffRate->bound[i][j][0] != -1)
		/* if no threshold or above threshold */
		{
		    offset += -(pRawDiffRate->bound)[i][j][0];
		    jLower = IMin(jLower, (pRawDiffRate->bound)[i][j][0]);
		    jUpper = IMax(jUpper, (pRawDiffRate->bound)[i][j][1]);
		    for (k = (pRawDiffRate->bound)[i][j][0]; 
			 k <= (pRawDiffRate->bound)[i][j][1]; k++)
		    {
#ifdef DEBUG
			CheckIndex(0, i+1, k, "FoldDiffRate");
			CheckIndex(0, pRawDiffRate->numberOfElements, k+offset,
				   "FoldDiffRate");
#endif
			(pDiffRate->diffRate)[i][k] += 
			    (pBgPhotonDensity->vector)[j]*
			    (pRawDiffRate->diffRate)[k+offset];
		    }
		    offset += (pRawDiffRate->bound)[i][j][1] + 1;
		    /* reset the offset index */
		}
	    }
	    (pDiffRate->bound)[i][0] = IMin((pDiffRate->bound)[i][0], jLower);
	    (pDiffRate->bound)[i][1] = IMax((pDiffRate->bound)[i][1], jUpper);
        }
    }
}


void FoldICS(const dCVector* pBgPhotonDensity, 
	     const RawTotalRate* ICSTotalRate,
	     const RawDiffRate* ICSPhotonRate, const RawDiffRate* ICSScatRate,
	     TotalRate* leptonTotalRate, DiffRate* leptonPhotonRate,
	     DiffRate* leptonScatRate)
{
    FoldTotalRate(pBgPhotonDensity, ICSTotalRate, leptonTotalRate);
    FoldDiffRate(pBgPhotonDensity, ICSScatRate, leptonScatRate, 1, 
        leptonTotalRate);
    FoldDiffRate(pBgPhotonDensity, ICSPhotonRate, leptonPhotonRate, 0);
}


void FoldTPP(const dCVector* pBgPhotonDensity, const dCVector* pEnergy,
	     const RawTotalRate* TPPTotalRate, const RawDiffRate* TPPDiffRate,
	     TotalRate* leptonTotalRate, DiffRate* leptonScatRate, 
	     DiffRate* leptonExchRate, dCVector* otherLoss)
{
    int i;
    int j;

    if ((pEnergy->dimension != TPPTotalRate->mainDimension) ||
	(TPPTotalRate->bgDimension != pBgPhotonDensity->dimension) ||
	(otherLoss->dimension != pEnergy->dimension))
    {
	Error("FoldTPP: inconsistent dimensions", PROGRAM_ERROR);
    }

    /* add TPP total rates in continuous energy loss */
    for (i = 0; i < TPPTotalRate->mainDimension; i++)
    {
        (otherLoss->vector)[i] = 0.;
        for (j = 0; j < TPPTotalRate->bgDimension; j++)
	{
            (otherLoss->vector)[i] += -(pEnergy->vector)[i]*
	        (pBgPhotonDensity->vector)[j]*
	        (TPPTotalRate->totalRate)[i][j];
	}
    }
    
    FoldDiffRate(pBgPhotonDensity, TPPDiffRate, leptonScatRate, 1,
        leptonTotalRate);
    FoldDiffRate(pBgPhotonDensity, TPPDiffRate, leptonExchRate, 0);
}


void FoldPP(const dCVector* pBgPhotonDensity, const RawTotalRate* PPTotalRate,
	    const RawDiffRate* PPDiffRate, TotalRate* photonTotalRate, 
	    DiffRate* photonLeptonRate)
{
    FoldTotalRate(pBgPhotonDensity, PPTotalRate, photonTotalRate);
    FoldDiffRate(pBgPhotonDensity, PPDiffRate, photonLeptonRate, 0);
}


void FoldDPP(const dCVector* pBgPhotonDensity, const RawTotalRate* DPPRate,
	     TotalRate* photonTotalRate, DiffRate* photonLeptonRate)
/* I adopt a simple model where one pair of e-/e+ gets all the energy
   equally (1/2); see PRD paper */
{
    int i;
    int j;
    int jLower;
    int jUpper;
    int num_main_bins;
    int num_bg_bins;

    double averaging_factor;
    int offset;
    double ratio;


    num_main_bins = DPPRate->mainDimension;
    num_bg_bins = DPPRate->bgDimension;

    if ((pBgPhotonDensity->dimension != num_bg_bins) ||
	(photonTotalRate->mainDimension != num_main_bins) ||
	(photonLeptonRate->mainDimension != num_main_bins))
    {
	Error("FoldDPP: inconsistent dimensions", PROGRAM_ERROR);
    }

    averaging_factor = (pow(10., 1./2./BINS_PER_DECADE) + 
	pow(10., -1./2./BINS_PER_DECADE))/2.;
    /* although this could be supplied through arguments, I provide a local
       version to keep the modality */
    offset = -(int)(BINS_PER_DECADE*log10(averaging_factor/2.) + 0.5);
    ratio = offset/BINS_PER_DECADE - log10(2.);

    for (i = 0; i < num_main_bins; i++)
    {
	jLower = photonLeptonRate->bound[i][0];
	jUpper = photonLeptonRate->bound[i][1];
        for (j = 0; j < num_bg_bins; j++)
        {
	    (photonTotalRate->totalRate)[i] += (pBgPhotonDensity->vector)[j]*
	        (DPPRate->totalRate)[i][j];
            if ((DPPRate->totalRate)[i][j] > 0.)
	    {
		/* this is supplanted by the new implementation...
                (photonLeptonRate->diffRate)[i][i-6] += 
                    (pBgPhotonDensity->vector)[j]*
		    (DPPRate->totalRate)[i][j]*(1. - alpha);
                (photonLeptonRate->diffRate)[i][i-7] += 
		    (pBgPhotonDensity->vector)[j]*
		    (DPPRate->totalRate)[i][j]*alpha;
		jLower = IMin((photonLeptonRate->bound)[i][0], i-7);
		jUpper = IMax((photonLeptonRate->bound)[i][1], i-6);
		*/
		if (ratio < 1.)
		{
		    if (i-offset >= 0)
		    {
#ifdef DEBUG
			CheckIndex(0, num_main_bins, i-offset, "FoldDPP");
#endif
			(photonLeptonRate->diffRate)[i][i-offset] +=
			    (pBgPhotonDensity->vector)[j]*
			    (DPPRate->totalRate)[i][j]*ratio;
			
			jLower = IMin(jLower, i-offset);
			jUpper = IMax(jUpper, i-offset);
		    }
		    if (i-offset-1 >= 0)
		    {
#ifdef DEBUG
			CheckIndex(0, num_main_bins, i-offset-1, "FoldDPP");
#endif
			(photonLeptonRate->diffRate)[i][i-offset-1] +=
			    (pBgPhotonDensity->vector)[j]*
			    (DPPRate->totalRate)[i][j]*(1. - ratio);
			
			jLower = IMin(jLower, i-offset-1);
		    }
		}
		else
		{
		    if (i-offset >= 0)
		    {
#ifdef DEBUG
			CheckIndex(0, num_main_bins, i-offset, "FoldDPP");
#endif
			(photonLeptonRate->diffRate)[i][i-offset] +=
		            (pBgPhotonDensity->vector)[j]*
			    (DPPRate->totalRate)[i][j]*(2. - ratio);

			jLower = IMin(jLower, i-offset);
			jUpper = IMax(jUpper, i-offset);
		    }
		    if (i-offset+1 >= 0)
		    {
#ifdef DEBUG
			CheckIndex(0, num_main_bins, i-offset+1, "FoldDPP");
#endif
			(photonLeptonRate->diffRate)[i][i-offset+1] +=
			    (pBgPhotonDensity->vector)[j]*
			    (DPPRate->totalRate)[i][j]*(ratio - 1.);

			jUpper = IMax(jUpper, i-offset+1);
		    }
		}
            }
        }
	(photonLeptonRate->bound)[i][0] = jLower;
	(photonLeptonRate->bound)[i][1] = jUpper;
	/* update bounds */
    }
}


void FoldPPPNucleon(const dCVector* pBgPhotonDensity, 
		    const RawTotalRate* PPPProtonLossRate,
		    const RawTotalRate* PPPNeutronLossRate, 
		    const RawDiffRate* PPPProtonScatRate,
		    const RawDiffRate* PPPProtonNeutronRate,
		    const RawDiffRate* PPPNeutronProtonRate,
		    TotalRate* protonTotalRate, TotalRate* neutronTotalRate,
		    DiffRate* protonScatRate, DiffRate* protonNeutronRate,
		    DiffRate* neutronProtonRate) 
{
    FoldTotalRate(pBgPhotonDensity, PPPProtonLossRate, protonTotalRate);
    FoldTotalRate(pBgPhotonDensity, PPPNeutronLossRate, neutronTotalRate);

    /*---- nucleon -> nucleon from PPP ----*/
    FoldDiffRate(pBgPhotonDensity, PPPProtonScatRate, protonScatRate, 
        2, protonTotalRate, neutronTotalRate);
    FoldDiffRate(pBgPhotonDensity, PPPProtonNeutronRate, protonNeutronRate, 
        0);
    FoldDiffRate(pBgPhotonDensity, PPPNeutronProtonRate, neutronProtonRate, 
        0);
}


void FoldPPPSecondary(const dCVector* pBgPhotonDensity,
		      const RawDiffRate* PPPProtonPhotonRate,
		      const RawDiffRate* PPPProtonElectronRate,
		      const RawDiffRate* PPPProtonPositronRate,
		      const RawDiffRate* PPPNeutronElectronRate,
		      const RawDiffRate* PPPProtonElectronNeutrinoRate,
		      const RawDiffRate* PPPProtonAntiElectronNeutrinoRate,
		      const RawDiffRate* PPPProtonMuonNeutrinoRate,
		      const RawDiffRate* PPPProtonAntiMuonNeutrinoRate,
		      const RawDiffRate* PPPNeutronAntiElectronNeutrinoRate,
		      const RawDiffRate* PPPNeutronMuonNeutrinoRate,
		      const RawDiffRate* PPPNeutronAntiMuonNeutrinoRate,
		      DiffRate* protonPhotonRate,
		      DiffRate* protonElectronRate, 
		      DiffRate* protonPositronRate,
		      DiffRate* neutronElectronRate, 
		      DiffRate* neutronPositronRate,
		      DiffRate* protonElectronNeutrinoRate,
		      DiffRate* protonAntiElectronNeutrinoRate,
		      DiffRate* protonMuonNeutrinoRate,
		      DiffRate* protonAntiMuonNeutrinoRate,
		      DiffRate* neutronAntiElectronNeutrinoRate,
		      DiffRate* neutronMuonNeutrinoRate,
		      DiffRate* neutronAntiMuonNeutrinoRate)
{
    /*---- nucleon -> EM species (gamma, e+/-) from PPP ----*/
    FoldDiffRate(pBgPhotonDensity, PPPProtonPhotonRate, protonPhotonRate, 
        0);
    FoldDiffRate(pBgPhotonDensity, PPPProtonPositronRate, protonPositronRate, 
        0);
    FoldDiffRate(pBgPhotonDensity, PPPNeutronElectronRate, 
        neutronElectronRate, 0);
    FoldDiffRate(pBgPhotonDensity, PPPProtonElectronRate, protonElectronRate, 
        0);
    FoldDiffRate(pBgPhotonDensity, PPPProtonElectronRate, neutronPositronRate,
        0);

    /*---- nucleons -> neutrinos from PPP ----*/
    FoldDiffRate(pBgPhotonDensity, PPPProtonMuonNeutrinoRate,
        protonMuonNeutrinoRate, 0);
    FoldDiffRate(pBgPhotonDensity, PPPNeutronAntiMuonNeutrinoRate,
        neutronAntiMuonNeutrinoRate, 0);
    FoldDiffRate(pBgPhotonDensity, PPPProtonAntiMuonNeutrinoRate,
        protonAntiMuonNeutrinoRate, 0);
    FoldDiffRate(pBgPhotonDensity, PPPNeutronMuonNeutrinoRate,
        neutronMuonNeutrinoRate, 0);
    FoldDiffRate(pBgPhotonDensity, PPPProtonElectronNeutrinoRate,
        protonElectronNeutrinoRate, 0);
    FoldDiffRate(pBgPhotonDensity, PPPNeutronAntiElectronNeutrinoRate,
        neutronAntiElectronNeutrinoRate, 0);
    FoldDiffRate(pBgPhotonDensity, PPPProtonAntiElectronNeutrinoRate,
        protonAntiElectronNeutrinoRate, 0);
}

void FoldNPPNucleon(const dCVector* pBgPhotonDensity, const dCVector* pEnergy,
		    const RawTotalRate* NPPTotalRate, 
		    dCVector* protonContinuousLoss)
{
    int i;
    int j;
    int num_main_bins;
    int num_bg_bins;

    num_main_bins = pEnergy->dimension;
    num_bg_bins = pBgPhotonDensity->dimension;

    if ((NPPTotalRate->mainDimension != num_main_bins) ||
	(NPPTotalRate->bgDimension != num_bg_bins) ||
	(protonContinuousLoss->dimension != num_main_bins))
    {
	Error("FoldNPP: inconsistent dimensions", PROGRAM_ERROR);
    }

    /*---- continuous energy loss ----*/
    for (i = 0; i < num_main_bins; i++)
    {
	(protonContinuousLoss->vector)[i] = 0.;
	/* this is very important! */
	for (j = 0; j < num_bg_bins; j++)
	{
	    (protonContinuousLoss->vector)[i] += -(pEnergy->vector)[i]*
	        (pBgPhotonDensity->vector)[j]*
	        (NPPTotalRate->totalRate)[i][j];
	}
    }
}

void FoldNPPSecondary(const dCVector* pBgPhotonDensity,
		      const RawDiffRate* NPPDiffRate,
		      DiffRate* protonPositronRate, 
		      DiffRate* protonElectronRate)
{
    FoldDiffRate(pBgPhotonDensity, NPPDiffRate, protonPositronRate, 0);
    FoldDiffRate(pBgPhotonDensity, NPPDiffRate, protonElectronRate, 0);
}


void MapNeutTotalRate(const double redshift, const int lastIndex,
		      const int tauNeutrinoMassSwitch, 
		      const TotalRate* NNTotalRate, TotalRate* pRate)
/* This is for m_nu = 0 */
{
    double redshiftFactor;
    int i;
    int offset;
    int iNew;
    int num_main_bins;

    num_main_bins = pRate->mainDimension;

    if (NNTotalRate->mainDimension != pRate->mainDimension)
    {
	Error("ManNeutTotalRate: inconsistent dimensions", PROGRAM_ERROR);
    }
    redshiftFactor = (1. + redshift)*(1. + redshift)*(1. + redshift);

    if (tauNeutrinoMassSwitch == 2)    /* massless case */
    {
	if (lastIndex == 0)    /* away from z = 0 (simple sliding applies) */
	{
	    offset = (int)(BINS_PER_DECADE*log10(1. + redshift));
	    for (i = 0; i < num_main_bins; i++)
	    {
		iNew = i + offset;
		/* map to the right index */
		if (iNew >= num_main_bins)
		    iNew = num_main_bins - 1;
		/* if the index is beyond range, simply set it to maximum
		   (unsatisfactory?) */
#ifdef DEBUG
		CheckIndex(0, num_main_bins, iNew, "MapNeutTotalRate");
#endif
		(pRate->totalRate)[i] = redshiftFactor*
		    (NNTotalRate->totalRate)[iNew];
	    }
	}
	else
	{
	    double fraction;

	    fraction = BINS_PER_DECADE*log10(1. + redshift);
	    offset = (int)fraction;
	    for (i = 0; i < num_main_bins; i++)
	    {
		iNew = i + offset;
		if (iNew < num_main_bins - 1)
		{
#ifdef DEBUG
		    CheckIndex(0, num_main_bins, iNew, "MapNeutTotalRate");
#endif
		    (pRate->totalRate)[i] = redshiftFactor*
		        ((NNTotalRate->totalRate)[iNew]*
			(1. - fraction) + (NNTotalRate->totalRate)[iNew+1]*
			 fraction);
		}
		else
		{
		    (pRate->totalRate)[i] = redshiftFactor*
		        (NNTotalRate->totalRate)[num_main_bins-1];
		}
	    }
	}
    }
    else    /* massive neutrinos */
    {
	for (i = 0; i < num_main_bins; i++)
	{
	    (pRate->totalRate)[i] = redshiftFactor*
	        (NNTotalRate->totalRate)[i];
	}
    }
}

void MapNeutDiffRate(const double redshift, const int lastIndex,
		     const int tauNeutrinoMassSwitch, 
		     const DiffRate* NNDiffRate, DiffRate* pRate)
{
    int i;
    int j;
    int offset;
    int iNew;
    int jNew;
    double redshiftFactor;
    int num_main_bins;

    num_main_bins = pRate->mainDimension;

    if (NNDiffRate->mainDimension != pRate->mainDimension)
    {
	Error("ManNeutDiffRate: inconsistent dimensions", PROGRAM_ERROR);
    }

    redshiftFactor = (1. + redshift)*(1. + redshift)*(1. + redshift);

    if (tauNeutrinoMassSwitch == 2)    /* massless case */
    {
	if (lastIndex == 0)    /* away from z = 0 (simple sliding applies) */
	{
	    offset = (int)(BINS_PER_DECADE*log10(1. + redshift));
	    for (i = 0; i < num_main_bins; i++)
	    {
		iNew = i + offset;
		/* map to the right index */
		if (iNew >= num_main_bins)
		    iNew = num_main_bins - 1;
		/* if the index is beyond range, simply set it to maximum
		   (unsatisfactory?) */
		for (j = 0; j < num_main_bins; j++)
		{
		    jNew = j + offset;
		    if (jNew >= num_main_bins)
			jNew = num_main_bins - 1;

#ifdef DEBUG
		    CheckIndex(0, num_main_bins, iNew, "MapNeutDiffRate");
		    CheckIndex(0, num_main_bins, jNew, "MapNeutDiffRate");
#endif
		    (pRate->diffRate)[i][j] = redshiftFactor*
		        (NNDiffRate->diffRate)[iNew][jNew];
		}
	    }
	}
	else
	{
	    double fraction;

	    fraction = BINS_PER_DECADE*log10(1. + redshift);
	    offset = (int)fraction;
	    for (i = 0; i < num_main_bins; i++)
	    {
		iNew = i + offset;
		for (j = 0; j < num_main_bins; j++)
		{
		    jNew = j + offset;
		    if (iNew < num_main_bins - 1)
		    {
#ifdef DEBUG
			CheckIndex(0, num_main_bins, iNew, "MapNeutDiffRate");
			CheckIndex(0, num_main_bins, jNew, "MapNeutDiffRate");
			CheckIndex(0, num_main_bins, iNew+1, 
				   "MapNeutDiffRate");
			CheckIndex(0, num_main_bins, jNew+1, 
				   "MapNeutDiffRate");
#endif
			(pRate->diffRate)[i][j] = redshiftFactor*
			    ((NNDiffRate->diffRate)[iNew][jNew]*
			    (1. - fraction)*(1. - fraction) + 
			    ((NNDiffRate->diffRate)[iNew+1][jNew] +
			    (NNDiffRate->diffRate)[iNew][jNew+1])*fraction*
			    (1. - fraction) +
			    (NNDiffRate->diffRate)[iNew+1][jNew+1]*fraction*
			    fraction);
		    }
		    else if (jNew < num_main_bins - 1)
		    {
#ifdef DEBUG
			CheckIndex(0, num_main_bins, jNew, "MapNeutDiffRate");
			CheckIndex(0, num_main_bins, jNew+1, 
				   "MapNeutDiffRate");
#endif
			(pRate->diffRate)[i][j] = redshiftFactor*
			    ((NNDiffRate->diffRate)[num_main_bins-1][jNew]*
			    (1. - fraction) + 
			    (NNDiffRate->diffRate)[num_main_bins-1][jNew+1]*
			    fraction);
		    }
		    else
		    {
			(pRate->diffRate)[i][j] = redshiftFactor*
			    (NNDiffRate->diffRate)[num_main_bins-1][num_main_bins-1];
		    }
		}
	    }
	}
    }
    else    /* massive neutrino: simply multiply by redshift factor */
    {
	for (i = 0; i < num_main_bins; i++)
	{
	    for (j = 0; j < num_main_bins; j++)
	    {
		(pRate->diffRate)[i][j] = redshiftFactor*
		    (NNDiffRate->diffRate)[i][j];
	    }
	}
    }

    /* take care of bounds: simply make it standard (0 <= k <= i) */
    for (i = 0; i < pRate->mainDimension; i++)
    {
	pRate->bound[i][0] = 0;
	pRate->bound[i][1] = i;
    }
}

void MapNeutRates(const double redshift, const int lastIndex,
		  const int tauNeutrinoMassSwitch,
		  const TotalRate* NNElNeutTotalRate,
		  const TotalRate* NNMuonNeutTotalRate,
		  const TotalRate* NNTauNeutTotalRate,
		  const DiffRate* NNElNeutScatRate, 
		  const DiffRate* NNElNeutMuonNeutRate,
		  const DiffRate* NNElNeutTauNeutRate,
		  const DiffRate* NNElNeutElectronRate,
		  const DiffRate* NNElNeutPhotonRate,
		  const DiffRate* NNElNeutProtonRate,
		  const DiffRate* NNMuonNeutElNeutRate,
		  const DiffRate* NNMuonNeutScatRate,
		  const DiffRate* NNMuonNeutTauNeutRate,
		  const DiffRate* NNMuonNeutElectronRate,
		  const DiffRate* NNMuonNeutPhotonRate,
		  const DiffRate* NNMuonNeutProtonRate,
		  const DiffRate* NNTauNeutElNeutRate, 
		  const DiffRate* NNTauNeutMuonNeutRate, 
		  const DiffRate* NNTauNeutScatRate, 
		  const DiffRate* NNTauNeutElectronRate, 
		  const DiffRate* NNTauNeutPhotonRate, 
		  const DiffRate* NNTauNeutProtonRate,
		  TotalRate* elNeutTotalRate, TotalRate* muonNeutTotalRate,
		  TotalRate* tauNeutTotalRate, DiffRate* elNeutScatRate, 
		  DiffRate* elNeutMuonNeutRate, DiffRate* elNeutTauNeutRate,
		  DiffRate* elNeutElectronRate, DiffRate* elNeutPhotonRate,
		  DiffRate* elNeutProtonRate, DiffRate* muonNeutElNeutRate,
		  DiffRate* muonNeutScatRate, DiffRate* muonNeutTauNeutRate,
		  DiffRate* muonNeutElectronRate, DiffRate* muonNeutPhotonRate,
		  DiffRate* muonNeutProtonRate, DiffRate* tauNeutElNeutRate, 
		  DiffRate* tauNeutMuonNeutRate, DiffRate* tauNeutScatRate, 
		  DiffRate* tauNeutElectronRate, DiffRate* tauNeutPhotonRate, 
		  DiffRate* tauNeutProtonRate)
{
    MapNeutTotalRate(redshift, lastIndex, tauNeutrinoMassSwitch, 
        NNElNeutTotalRate, elNeutTotalRate);
    MapNeutTotalRate(redshift, lastIndex, tauNeutrinoMassSwitch, 
        NNMuonNeutTotalRate, muonNeutTotalRate);
    MapNeutTotalRate(redshift, lastIndex, tauNeutrinoMassSwitch,
        NNTauNeutTotalRate, tauNeutTotalRate);
    MapNeutDiffRate(redshift, lastIndex, tauNeutrinoMassSwitch,
        NNElNeutScatRate, elNeutScatRate);
    MapNeutDiffRate(redshift, lastIndex, tauNeutrinoMassSwitch,
        NNElNeutMuonNeutRate, elNeutMuonNeutRate);
    MapNeutDiffRate(redshift, lastIndex, tauNeutrinoMassSwitch,
        NNElNeutTauNeutRate, elNeutTauNeutRate);
    MapNeutDiffRate(redshift, lastIndex, tauNeutrinoMassSwitch,
        NNElNeutElectronRate, elNeutElectronRate);
    MapNeutDiffRate(redshift, lastIndex, tauNeutrinoMassSwitch,
        NNElNeutPhotonRate, elNeutPhotonRate);
    MapNeutDiffRate(redshift, lastIndex, tauNeutrinoMassSwitch,
        NNElNeutProtonRate, elNeutProtonRate);
    MapNeutDiffRate(redshift, lastIndex, tauNeutrinoMassSwitch,
        NNMuonNeutElNeutRate, muonNeutElNeutRate);
    MapNeutDiffRate(redshift, lastIndex, tauNeutrinoMassSwitch,
        NNMuonNeutScatRate, muonNeutScatRate);
    MapNeutDiffRate(redshift, lastIndex, tauNeutrinoMassSwitch,
        NNMuonNeutTauNeutRate, muonNeutTauNeutRate);
    MapNeutDiffRate(redshift, lastIndex, tauNeutrinoMassSwitch,
        NNMuonNeutElectronRate, muonNeutElectronRate);
    MapNeutDiffRate(redshift, lastIndex, tauNeutrinoMassSwitch,
        NNMuonNeutPhotonRate, muonNeutPhotonRate);
    MapNeutDiffRate(redshift, lastIndex, tauNeutrinoMassSwitch,
        NNMuonNeutProtonRate, muonNeutProtonRate);
    MapNeutDiffRate(redshift, lastIndex, tauNeutrinoMassSwitch,
        NNTauNeutElNeutRate, tauNeutElNeutRate);
    MapNeutDiffRate(redshift, lastIndex, tauNeutrinoMassSwitch,
        NNTauNeutMuonNeutRate, tauNeutMuonNeutRate);
    MapNeutDiffRate(redshift, lastIndex, tauNeutrinoMassSwitch,
        NNTauNeutScatRate, tauNeutScatRate);
    MapNeutDiffRate(redshift, lastIndex, tauNeutrinoMassSwitch,
        NNTauNeutElectronRate, tauNeutElectronRate);
    MapNeutDiffRate(redshift, lastIndex, tauNeutrinoMassSwitch,
        NNTauNeutPhotonRate, tauNeutPhotonRate);
    MapNeutDiffRate(redshift, lastIndex, tauNeutrinoMassSwitch,
        NNTauNeutProtonRate, tauNeutProtonRate);
}

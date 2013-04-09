
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "dint/const.h"
#include "dint/cvector.h"
#include "dint/utilities.h"
#include "dint/spectrum.h"
#include "dint/frag.h"
#include "dint/decay.h"


void SetEnergyBins(const int min_energy_exp, dCVector* pEnergy, 
		   dCVector* pEnergyWidth)
{
    int i;
    double exponent;
    double temp;
    int num_bins;
    double averaging_factor;
    double binning_factor;

    num_bins = pEnergy->dimension;
    averaging_factor = (pow(10., 1./2./BINS_PER_DECADE) +
        pow(10., -1./2./BINS_PER_DECADE))/2.;
    binning_factor = pow(10., 1./2./BINS_PER_DECADE) -
        pow(10., -1./2./BINS_PER_DECADE);

    if (num_bins != pEnergyWidth->dimension)
    {
	Error("SetEnergyBins: inconsistent dimensions", PROGRAM_ERROR);
    }

    for (i = 0; i < num_bins; i++)
    {
        exponent = (double)min_energy_exp + (double)i/BINS_PER_DECADE;
        temp = pow(10., exponent)/ELECTRON_MASS;
        (pEnergy->vector)[i] = temp*averaging_factor;
        (pEnergyWidth->vector)[i] = temp*binning_factor;
    }
}

void SetDeltaG(const dCVector* pEnergy, dCVector* pDeltaG)
{
    int i;
    int num_main_bins;

    num_main_bins = pEnergy->dimension;

    if (pDeltaG->dimension != num_main_bins)
    {
	Error("SetDeltaG: inconsistent dimensions", PROGRAM_ERROR);
    }

    for (i = 0; i < num_main_bins - 1; i++)
    {
        (pDeltaG->vector)[i] = (pEnergy->vector)[i+1] - (pEnergy->vector)[i];
    }
    (pDeltaG->vector)[num_main_bins-1] = (pEnergy->vector)[num_main_bins-1]*
        pow(10., 1./BINS_PER_DECADE);
}

/* NOTE: this function is not used any more... */
void GetModeOfInput(FILE* input, int* pInputMethodSwitch)
{
    /* give a choice between manually entering the data and reading it
        from a data file */
    printf("        -----------------< Enter Model Parameters >");
    printf("-----------------\n\n");
    printf("1. Read from a file (1)/Manually enter parameters (0): ");
    scanf("%i", pInputMethodSwitch);

    if (*pInputMethodSwitch == 1)
    /* From data file */
    {
        char fileName[20];

        /* get the name of the file */
        printf("2. Name of parameter file: ");
	scanf("%s", fileName);
        input = SafeFOpen(fileName, "r");
    }
}

void BasicParameterInput(FILE* input, const int argc, int* pMinEnergyExp,
			 int* pNumSmallSteps, 
			 double* pConvergeParameter,
			 double* pStartingRedshift, 
			 double* pStartingDistanceInMpc)
{
    if (argc == 3)
    /* From data file */
    {
	printf("\n[Parameters read from the file]\n");
        /* read in parameters */
	fscanf(input, "%i", pMinEnergyExp);
	printf("1. Minimum energy: 10^%i eV\n", *pMinEnergyExp);
        fscanf(input, "%i", pNumSmallSteps);
	printf("2. Number of small steps: %i\n", *pNumSmallSteps);
        fscanf(input, "%lf", pConvergeParameter);
	printf("3. Delta (convergence parameter): %12.3E\n", 
            *pConvergeParameter);
        fscanf(input, "%lf", pStartingRedshift);
	if (*pStartingRedshift != 0.)
	{
	    printf("4. Maximal redshift: %10.3f\n", *pStartingRedshift);
	}
        fscanf(input, "%lf", pStartingDistanceInMpc);
	if (*pStartingDistanceInMpc >= 2.*C/3./H_0*1.e-5)
	{
	    Error("BasicParameterInput: distance larger than horizon",
		  PROGRAM_ERROR);
	}
	if (*pStartingDistanceInMpc != 0.)
	{
	    printf("4. Maximal distance (Mpc): %10.3f\n", 
                *pStartingDistanceInMpc);
	}
    }
    else if (argc == 2)
    /* Enter manually */
    {
	printf("1. Minimum energy exponent (e.g. 10^8 eV: 8): ");
	scanf("%i", pMinEnergyExp);
        printf("2. Number of small steps (suggestion: 6): ");
        scanf("%i", pNumSmallSteps);
        printf("3. Delta (convergence parameter) (suggestion: 1.e-6): ");
        scanf("%lf", pConvergeParameter);
        printf("4. Maximal redshift ");
	printf("(enter 0 if you will type in distance): ");
        scanf("%lf", pStartingRedshift);
        printf("5. Maximal distance in Mpc ");
	printf("(enter 0 if you typed in redshift): ");
        scanf("%lf", pStartingDistanceInMpc);
	if (*pStartingDistanceInMpc >= 2.*C/3./H_0*1.e-5)
	{
	    Error("BasicParameterInput: distance larger than horizon",
		  PROGRAM_ERROR);
	}
    }

    if (*pStartingDistanceInMpc == 0.)
    {
        *pStartingDistanceInMpc = 2.*C/3./H_0/1.e5*(1. - 
            pow(1. + *pStartingRedshift, -3./2.));
	printf("\n");
	printf("Maximal distance (Mpc): %12.6f\n", *pStartingDistanceInMpc);
    }
    if (*pStartingRedshift == 0.)
    {
        *pStartingRedshift = 
            pow(1. - 3.*H_0*1.e5*(*pStartingDistanceInMpc)/2./C, -2./3.) - 1.;
	printf("\n");
	printf("Maximal redshift: %12.6f\n", *pStartingRedshift);
    }
}

void InteractionParameterInput(FILE* input, const int argc,
			       int* pSynchrotronSwitch, double* pB_0, 
			       int* pTauNeutrinoMassSwitch, int* pICSSwitch, 
			       int* pPPSwitch, int* pTPPSwitch, 
			       int* pDPPSwitch, int* pPPPSwitch, 
			       int* pNPPSwitch, int* pNeutronDecaySwitch,
			       int* pNucleonToSecondarySwitch,
			       int* pNeutrinoNeutrinoSwitch)
{
    if (argc == 3)
    /* From data file */
    {
        fscanf(input, "%i", pSynchrotronSwitch);
	fscanf(input, "%lf", pB_0);
        if (*pSynchrotronSwitch == 1)
        {
	    printf("5. Synchrotron loss on -> B_0: %12.3E\n", *pB_0);
        }
        else
        {
            *pB_0 = 0;
	    printf("5. Synchrotron loss off\n");
        }
	fscanf(input, "%i", pTauNeutrinoMassSwitch);
	if (*pTauNeutrinoMassSwitch == 0)
	    printf("6. m_tau = m_e = m_nu = 1 eV\n");
	else
	    printf("6. m_tau = m_nu = 1 eV, m_e = 0.1 eV\n");

	/* read in interaction switches */
	printf("\n7. ");
	fscanf(input, "%i", pICSSwitch);
	if (*pICSSwitch != 0)
	    printf("ICS on... ");
	fscanf(input, "%i", pPPSwitch);
	if (*pPPSwitch != 0)
	    printf("PP on... ");
	fscanf(input, "%i", pTPPSwitch);
	if (*pTPPSwitch != 0)
	    printf("TPP on... ");
	fscanf(input, "%i", pDPPSwitch);
	if (*pDPPSwitch != 0)
	    printf("DPP on... ");
	fscanf(input, "%i", pPPPSwitch);
	if (*pPPPSwitch != 0)
	    printf("Photopion production on... ");
	fscanf(input, "%i", pNPPSwitch);
	if (*pNPPSwitch != 0)
	    printf("Proton pair production on... ");
	fscanf(input, "%i", pNeutronDecaySwitch);
	if (*pNeutronDecaySwitch != 0)
	    printf("Neutron decay on... ");
	fscanf(input, "%i", pNucleonToSecondarySwitch);
	if (*pNucleonToSecondarySwitch != 0)
	    printf("Nucleon secondary tables included... ");
	fscanf(input, "%i", pNeutrinoNeutrinoSwitch);
	if (*pNeutrinoNeutrinoSwitch != 0)
	    printf("Neutrino-neutrino interaction on... ");
	printf("\n");
    }
    else if (argc == 2)
    /* Enter manually */
    {
        printf("6. Synchrotron on/off? (on = 1, off = 0) ");
        scanf("%i", pSynchrotronSwitch);
        if (*pSynchrotronSwitch == 1)
        {
            printf("  6-1. Magnetic field (G): ");
            scanf("%lf", pB_0);
        }
        else
        {
            *pB_0 = 0;
        }

	printf("7. Neutrino mass (0: equally massive (1 eV), 1: tau = nu = 1, e = 0.1 eV): ");
	scanf("%lf", pTauNeutrinoMassSwitch);

	//	printf("\n\nNow turning to individual interactions...\n");
	printf("\n");
	printf("8-1. Inverse Compton scattering (ICS) on? (on = 1, off = 0) ");
	scanf("%i", pICSSwitch);
	printf("8-2. Pair production (PP) on? ");
	scanf("%i", pPPSwitch);
	printf("8-3. Triplet pair production (TPP) on? ");
	scanf("%i", pTPPSwitch);
	printf("8-4. Double pair production (DPP) on? ");
	scanf("%i", pDPPSwitch);
	printf("8-5. Photopion production on? ");
	scanf("%i", pPPPSwitch);
	printf("8-6. Proton pair production on? ");
	scanf("%i", pNPPSwitch);
	printf("8-7. Neutron decay on? ");
	scanf("%i", pNeutronDecaySwitch);
	printf("8-8. Nucleon secondary tables included? ");
	scanf("%i", pNucleonToSecondarySwitch);
	printf("8-9. Neutrino-neutrino interaction on? ");
	scanf("%i", pNeutrinoNeutrinoSwitch);
    }

    /* quick check whether secondary table switch is consistent */
    if (*pNucleonToSecondarySwitch == 1)
    {
        if ((*pPPPSwitch == 0) && (*pNPPSwitch == 0) && 
	    (*pNeutronDecaySwitch == 0))
	{
	    printf("WARNING: secondary table on when all interactions off?\n");
	}
    }
}

void ModelParameterInput(FILE* input, const int argc,
			 int* pSourceTypeSwitch, double* pMinDistance,
			 double* pBrightPhaseExp, int* pModelTypeSwitch)
{
    if (argc == 3)
    /* From data file */
    {
        fscanf(input, "%i", pSourceTypeSwitch);
	fscanf(input, "%lf", pMinDistance);
	fscanf(input, "%lf", pBrightPhaseExp);
	if (*pSourceTypeSwitch == 1)
	  {
	    printf("8. Source type: diffuse\n");
	    printf("8-1. Minimal distance: %12.6f\n", *pMinDistance);
	    printf("8-2. Bright phase exponent: %12.6f\n", *pBrightPhaseExp);
	  }
	else
	    printf("8. Source type: single\n");

	fscanf(input, "%i", pModelTypeSwitch);
	if (*pModelTypeSwitch == 0)
	    printf("9. Injection model: photon monoenergetic injection\n");
	else if (*pModelTypeSwitch == 1)
	    printf("9. Injection model: electron monoenergetic injection\n");
	else if (*pModelTypeSwitch == 2)
	    printf("9. Injection model: positron monoenergetic injection\n");
    }
    else if (argc == 2)
    /* Enter manually */
    {
        printf("9. Single source or diffuse sources? ");
	printf("(single = 0, diffuse = 1) ");
        scanf("%i", pSourceTypeSwitch);
	*pMinDistance = 0.;
	*pBrightPhaseExp = 1.5;
	if (*pSourceTypeSwitch == 1)
	  {
	    printf("9-1. Minimal distance: ");
	    scanf("%lf", pMinDistance);
	    printf("9-2. Bright phase exponent: ");
	    scanf("%lf", pBrightPhaseExp);
	  }
	printf("10. Injection model\n");
	printf("(photon = 0, electron = 1,\n");
	printf(" positron = 2; monoenergetic injection ? ");
	scanf("%i", pModelTypeSwitch);
    }

    printf("\n\nParameter input complete.\n");
}

void PrepareSpectra(const int sourceTypeSwitch, const Spectrum* pQ_0, 
                    Spectrum* pSpectrum, Spectrum* pSpectrumNew, 
		    Spectrum* pDerivative)
{
    if (pSpectrumNew->numberOfMainBins != pDerivative->numberOfMainBins)
    {
	Error("PrepareSpectra: inconsistent dimensions", PROGRAM_ERROR);
    }

    //    InitializeSpectrum(spectrumOld);
    if (sourceTypeSwitch == 0)  /* single source */
        SetSpectrum(pSpectrum, pQ_0);
    else    /* diffuse sources */
        InitializeSpectrum(pSpectrum);

    SetSpectrum(pSpectrumNew, pSpectrum);
    InitializeSpectrum(pDerivative);
}


void ComputeTotalInitialContent(const dCVector* pEnergy, const Spectrum* pQ_0, 
                                double* initialPhotonEnergy,
                                double* initialLeptonEnergy, 
				double* initialNucleonEnergy,
				double* initialNeutrinoEnergy,
                                double* initialTotalEnergy,
                                double* initialPhotonNumber,
                                double* initialLeptonNumber,
				double* initialNucleonNumber,
				double* initialNeutrinoNumber,
                                double* initialTotalNumber)
{
    int i;

    *initialPhotonEnergy = 0;
    *initialLeptonEnergy = 0.;
    *initialNucleonEnergy = 0.;
    *initialNeutrinoEnergy = 0.;
    *initialTotalEnergy = 0.;
    *initialPhotonNumber = 0;
    *initialLeptonNumber = 0.;
    *initialNucleonNumber = 0.;
    *initialNeutrinoNumber = 0.;
    *initialTotalNumber = 0.;
    
    for (i = 0; i < pEnergy->dimension; i++)
    {
        *initialPhotonEnergy += (pQ_0->spectrum)[PHOTON][i]*
	    (pEnergy->vector)[i];
        *initialPhotonNumber += (pQ_0->spectrum)[PHOTON][i];
        *initialLeptonEnergy += ((pQ_0->spectrum)[ELECTRON][i] + 
	    (pQ_0->spectrum)[POSITRON][i])*(pEnergy->vector)[i];
        *initialLeptonNumber += (pQ_0->spectrum)[ELECTRON][i] + 
	    (pQ_0->spectrum)[POSITRON][i];
    }
    *initialNucleonEnergy += GetNucleonEnergy(pQ_0, pEnergy);
    *initialNucleonNumber += GetNucleonNumber(pQ_0);
    *initialNeutrinoEnergy += GetNeutrinoEnergy(pQ_0, pEnergy);
    *initialNeutrinoNumber += GetNeutrinoNumber(pQ_0);
    *initialTotalEnergy = *initialPhotonEnergy + *initialLeptonEnergy +
        *initialNucleonEnergy + *initialNeutrinoEnergy;
    *initialTotalNumber = *initialPhotonNumber + *initialLeptonNumber +
        *initialNucleonNumber + *initialNeutrinoNumber;
}


void ComputeContinuousEnergyLoss(const int synchrotronSwitch, 
                                 const dCVector* synchrotronLoss, 
                                 const dCVector* otherLoss, 
                                 dCVector* continuousLoss)
{
    int i;

    if ((synchrotronLoss->dimension != otherLoss->dimension) ||
	(otherLoss->dimension != continuousLoss->dimension))
    {
	Error("ComputeContinuousEnergyLoss: inconsistent dimensions",
	      PROGRAM_ERROR);
    }

    for (i = 0; i < continuousLoss->dimension; i++)
    {
        (continuousLoss->vector)[i] = (otherLoss->vector)[i];
        if (synchrotronSwitch == 1)     /* synchrotron on */
	{
            (continuousLoss->vector)[i] += (synchrotronLoss->vector)[i];
	}
    }
}

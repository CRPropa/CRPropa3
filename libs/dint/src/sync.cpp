
#include "dint/sync.h"

void LoadSyncTable(dCVector* syncTable, string aDirTables)
{
  FILE* syncData;
  double temporary;
  int i;
  
  syncData = SafeFOpen((aDirTables+"/syncTable.dat").c_str(), "r");
  for (i = 0; i < syncTable->dimension; i++) {
    fscanf(syncData, "%lf %lf", &temporary, &((syncTable->vector)[i]));
  }
  fclose(syncData) ; // added June 2005 - E.A.

}


void InitializeSynchrotron(const double B_0, const dCVector* pEnergy, 
			   const dCVector* pEnergyWidth,
			   dCVector* synchrotronLoss, DiffRate* syncRate,
			   string aDirTables)
/* this function calculates the synchrotron loss rate for the electron
    and the synchrotron radiation spectrum for photons
    The radom direction of the magnetic field requires me to do an
    average over the pitch angle, and it results in an extra factor of
    2/3, and slightly modifies the functional form of F(x) from the
    purely perpendicular case   */
/* The above comment is from the stand alone version. Compared to this version the constants seem to be fixed to correct for a perpendicular field as input. JK, 2011*/
{
    const int SYNC_TABLE_SIZE = 161;
    int i;
    double offset;    /* parameter to compute bounds correctly */
    dCVector syncTable;
    double x0;

    double ECrit;
    double EMin;
    double EMax;
    int jMin;
    int jMax;
    int j;
    double xx;
    double function;
    int iTable;
    double eLower;
    double eUpper;
    int num_main_bins;

    num_main_bins = pEnergy->dimension;

    if ((num_main_bins != pEnergyWidth->dimension) || 
	(num_main_bins != syncRate->mainDimension) ||
	(num_main_bins != synchrotronLoss->dimension))
    {
	Error("InitializeSynchrotron: inconsistent dimensions", PROGRAM_ERROR);
    }

    New_dCVector(&syncTable, SYNC_TABLE_SIZE);

    /* read in table */

    LoadSyncTable(&syncTable, aDirTables);

    /* calculate the rates */
    offset = BINS_PER_DECADE*(log10(ELECTRON_MASS) - MAX_ENERGY_EXP) +
        num_main_bins + 0.5;
    x0 = B_0*5.86667629e-4*DISTANCE_UNIT;
    /* x0: eB/m_e in inverse pc */

    for (i = 0; i < num_main_bins; i++)
    {
        /* electron loss rate */
        (synchrotronLoss->vector)[i] = -B_0*B_0*(pEnergy->vector)[i]*
	    (pEnergy->vector)[i]*2.*4.86037e4*VOLUME_UNIT;
        
        /* calculate the radiation spectrum */
        ECrit = 3./2.*B_0/4.414034e13*(pEnergy->vector)[i]*
	    (pEnergy->vector)[i];
        /* ECrit: critical energy in electron mass */

        EMin = 2./3./((pEnergy->vector)[i]*(pEnergy->vector)[i]*
            (pEnergy->vector)[i])*ECrit;
        EMax = 5.*ECrit;
        /* EMin/EMax: useful range for x (in electron mass) */

        /* set up the range for photons -> determine bounds */
	jMin = IMax((int)(BINS_PER_DECADE*log10(EMin) + offset), 0);
	jMax = IMin((int)(BINS_PER_DECADE*log10(EMax) + offset),
	    num_main_bins - 1);
        if (jMax >= 0)  /* normal case */
        {
            (syncRate->bound)[i][0] = jMin;
            (syncRate->bound)[i][1] = jMax;

            for (j = (syncRate->bound)[i][0]; 
		 j <= (syncRate->bound)[i][1]; j++)
            {
                xx = (pEnergy->vector)[j]/ECrit;

                if (log10(xx) < -7.)
                {
                    function = pow(xx, 1./3.)*2.1809736;
		    /* use an approximation for x below 10^-7 */
                }
                else
                {
		  /* compute F(x) at given photon energy by linear
		     extrapolation in logarithm */
                    iTable = (int)(log10(xx)*20. + 140);
                    eLower = (double)(iTable - 140)/20.;
                    eUpper = (double)(iTable - 139)/20.;
		    /* numbers 140, 139, and 20 are properties of table */
                    function = exp(log(syncTable.vector[iTable]) + (log10(xx) -
                        eLower)/(eUpper - eLower)*
                        (log(syncTable.vector[iTable+1]) - 
			 log(syncTable.vector[iTable])));
                }
                (syncRate->diffRate)[i][j] = 2.756644477e-1/137.036*x0*
		    function/(pEnergy->vector)[j]*(pEnergyWidth->vector)[j];
                /* 2.7566...: sqtr(3)/(2pi), 1/137: e^2, x0: eB/m_e */
            }
        }
    }

    Delete_dCVector(&syncTable) ;

}

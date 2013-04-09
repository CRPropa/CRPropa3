
#ifndef _SYNC_H_
#define _SYNC_H_

#include "dint/rate.h"
#include "dint/cvector.h"
#include <stdio.h>
#include <math.h>
#include "dint/utilities.h"
#include "dint/const.h"
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
using namespace std;

void LoadSyncTable(dCVector* syncTable, string aDirTables);
void InitializeSynchrotron(const double B_0, const dCVector* pEnergy, 
			   const dCVector* pEnergyWidth,
			   dCVector* synchrotronLoss, DiffRate* syncRate,
			   string aDirTables);

#endif

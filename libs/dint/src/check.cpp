#include <stdio.h>
#include <stdlib.h>
#include "dint/error.h"
#include "dint/cvector.h"
#include "dint/const.h"


void DumpArray(const dCVector* pVector)
{
    int i;
    for (i = 0; i < pVector->dimension; i++)
    {
	printf("%15.6E\n", pVector->vector[i]*ELECTRON_MASS);
    }
}

void CheckIndex(const int lowerLimit, const int upperLimit, const int i,
		const char* functionName)
{
    if (i < lowerLimit || i >= upperLimit)
    {
	printf("%s: index out of bounds!!!\n", functionName);
	printf("Should have satisfied %i <= %i < %i.\n", lowerLimit, i,
	       upperLimit);
	exit (ARRAY_ERROR);
    }
}

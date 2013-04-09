#include <stdio.h>
#include "dint/rate.h"
#include "dint/vector.h"
#include "dint/utilities.h"


/* this is actually the class constructor for RawTotalRate */
void NewRawTotalRate(RawTotalRate* pRate, const int num_main_bins,
		     const int num_bg_bins)
{
    pRate->mainDimension = num_main_bins;
    pRate->bgDimension = num_bg_bins;
    pRate->totalRate = New_dMatrix(num_main_bins, num_bg_bins);

    InitializeRawTotalRate(pRate);
}

void DeleteRawTotalRate(RawTotalRate* pRate)
{
    Delete_dMatrix(pRate->totalRate);
}

void InitializeRawTotalRate(RawTotalRate* pRate)
{
    int i;
    int j;

    for (i = 0; i < pRate->mainDimension; i++)
    {
	for (j = 0; j < pRate->bgDimension; j++)
	{
	    pRate->totalRate[i][j] = 0.;
	}
    }
}

void CopyRawTotalRate(RawTotalRate* pLRate, const RawTotalRate* pRRate)
{
    int i;
    int j;

    pLRate->mainDimension = pRRate->mainDimension;
    pLRate->bgDimension = pRRate->bgDimension;

    Delete_dMatrix(pLRate->totalRate);

    pLRate->totalRate = New_dMatrix(pRRate->mainDimension, 
	pRRate->bgDimension);

    for (i = 0; i < pRRate->mainDimension; i++)
    {
	for (j = 0; j < pRRate->bgDimension; j++)
	{
	    pLRate->totalRate[i][j] = pRRate->totalRate[i][j];
	}
    }
    /* deep copy */
}

void ClipRawTotalRate(RawTotalRate* pRate, const int newSize)
{
    int i;
    int j;
    int clip;
    RawTotalRate tempRate;

    printf("Clipping RawTotalRate.  From %5i...  ", pRate->mainDimension);
    /* let user know */

    clip = pRate->mainDimension - newSize;
    NewRawTotalRate(&tempRate, newSize, pRate->bgDimension);

    for (i = 0; i < newSize; i++)
    {
	int iOld = i + clip;
	for (j = 0; j < pRate->bgDimension; j++)
	{
#ifdef DEBUG
	    CheckIndex(0, pRate->mainDimension, iOld, "ClipRawTotalRate");
#endif
	    tempRate.totalRate[i][j] = pRate->totalRate[iOld][j];
	}
    }

    CopyRawTotalRate(pRate, &tempRate);
    DeleteRawTotalRate(&tempRate);
    
    printf("to %5i.\n", pRate->mainDimension);
}

void EnlargeRawTotalRate(RawTotalRate* pRate, const int newSize)
{
    int i;
    int j;
    RawTotalRate tempRate;
    int numAddedCells;


    numAddedCells = newSize - pRate->mainDimension;
    NewRawTotalRate(&tempRate, newSize, pRate->bgDimension);

    for (i = 0; i < newSize; i++)
    {
	int iOld = i - numAddedCells;
	for (j = 0; j < pRate->bgDimension; j++)
	{
	    if (iOld >= 0)    /* range where original matrix is */
	    {
#ifdef DEBUG
		CheckIndex(0, pRate->mainDimension, iOld, 
			   "EnlargeRawTotalRate");
#endif
		tempRate.totalRate[i][j] = pRate->totalRate[iOld][j];
	    }
	    else    /* enlarged area is filled with 0 for convenience */
	    {
		tempRate.totalRate[i][j] = 0.;
	    }
	}
    }
    
    CopyRawTotalRate(pRate, &tempRate);
    DeleteRawTotalRate(&tempRate);
}

void ModifyRawTotalRate(RawTotalRate* pRate, const int newSize)
{
    if (pRate->mainDimension > newSize)
    {
	ClipRawTotalRate(pRate, newSize);
    }
    else if (pRate->mainDimension < newSize)
    {
	EnlargeRawTotalRate(pRate, newSize);
    }
    else
    {
	/* do nothing */
    }
}

void ReadRawTotalRate(RawTotalRate* pRate, const char* filename)
{
    FILE* file;

    file = SafeFOpen(filename, "r");
    binfread((pRate->totalRate)[0], sizeof(double), (pRate->mainDimension)*
        (pRate->bgDimension), file);
    fclose(file);
}



/* this is actually the class constructor for RawDiffRate */
void NewRawDiffRate(RawDiffRate* pRate, const int num_main_bins, 
		    const int num_bg_bins, const int num_elements)
{
    pRate->mainDimension = num_main_bins;
    pRate->bgDimension = num_bg_bins;
    pRate->numberOfElements = num_elements;
    pRate->bound = New_i3Tensor(num_main_bins, num_bg_bins, 2);
    pRate->diffRate = New_dVector(num_elements);

    InitializeRawDiffRate(pRate);
}

void DeleteRawDiffRate(RawDiffRate* pRate)
{
    Delete_dVector(pRate->diffRate);
    Delete_i3Tensor(pRate->bound);
}

void InitializeRawDiffRate(RawDiffRate* pRate)
{
    int i;
    int j;

    for (i = 0; i < pRate->mainDimension; i++)
    {
	for (j = 0; j < pRate->bgDimension; j++)
	{
	    pRate->bound[i][j][0] = -1;
	    pRate->bound[i][j][1] = -1;
	}
    }
    /* note how the bounds are invalidated */

    for (i = 0; i < pRate->numberOfElements; i++)
    {
	pRate->diffRate[i] = 0.;
    }
}

void CheckRawDiffRate(RawDiffRate* pRate)
{
    int i;
    int j;
    int counter;

    for (counter = 0, i = 0; i < pRate->mainDimension; i++)
    {
	for (j = 0; j < pRate->bgDimension; j++)
	{
	    if (pRate->bound[i][j][0] != -1)
	    {
		counter += pRate->bound[i][j][1] - pRate->bound[i][j][0] + 1;
	    }
	}
    }

    if (counter != pRate->numberOfElements)
    {
	Error("CheckRawDiffRate: inconsistent rate", PROGRAM_ERROR);
    }
}

#ifdef DEBUG
double RawDiffRateElement(const RawDiffRate* pRate, const int i, const int j,
			  const int k)
{
    int n;
    int p;
    int offset;


    if ((i < 0) || (j < 0) || (i >= pRate->mainDimension) ||
	(j >= pRate->bgDimension))
    {
	printf("i: %3i, j: %3i\n", i, j);
	Error("RawDiffRateElement: i or j out of bounds", PROGRAM_ERROR);
    }
    if ((pRate->bound)[i][j][0] == -1)    /* out of range */
    {
	printf("i: %3i, j: %3i\n", i, j);
	Error("RawDiffRateElement: i or j out of range", PROGRAM_ERROR);
    }
    else
    {
	if ((k < (pRate->bound)[i][j][0]) || (k > (pRate->bound)[i][j][1]))
	{
	    printf("i, j, k: %3i %3i %3i\n", k);
	    Error("RawDiffRateElement: k out of range", PROGRAM_ERROR);
	}
    }
    
    offset = 0;
    for (n = 0; n < i; n++)
    {
	for (p = 0; p < pRate->bgDimension; p++)
	{
	    if ((pRate->bound)[n][p][0] != -1)
	    {
		offset += ((pRate->bound)[n][p][1] - (pRate->bound)[n][p][0] +
		    1);
	    }
	}
    }
    for (p = 0; p < j; p++)
    {
	offset += ((pRate->bound)[i][p][1] - (pRate->bound)[i][p][0] + 1);
    }
    offset += -(pRate->bound)[i][j][0] + k;
    if (offset >= pRate->numberOfElements || offset < 0)
    {
	printf("offset: %7i, number of elements: %7i\n", offset,
	       pRate->numberOfElements);
	Error("RawDiffRateElement: main index out of range", PROGRAM_ERROR);
    }

    return (pRate->diffRate)[offset];
}
#endif

void CopyRawDiffRate(RawDiffRate* pLRate, const RawDiffRate* pRRate)
{
    int i;
    int j;

    pLRate->mainDimension = pRRate->mainDimension;
    pLRate->bgDimension = pRRate->bgDimension;
    pLRate->numberOfElements = pRRate->numberOfElements;

    Delete_dVector(pLRate->diffRate);
    Delete_i3Tensor(pLRate->bound);

    pLRate->bound = New_i3Tensor(pRRate->mainDimension, pRRate->bgDimension,
        2);
    pLRate->diffRate = New_dVector(pRRate->numberOfElements);

    for (i = 0; i < pRRate->mainDimension; i++)
    {
	for (j = 0; j < pRRate->bgDimension; j++)
	{
	    pLRate->bound[i][j][0] = pRRate->bound[i][j][0];
	    pLRate->bound[i][j][1] = pRRate->bound[i][j][1];
	}
    }
    for (i = 0; i < pRRate->numberOfElements; i++)
    {
	pLRate->diffRate[i] = pRRate->diffRate[i];
    }
}

void CopyRawDiffRateBound(RawDiffRate* pLRate, const RawDiffRate* pRRate)
{
    int i;
    int j;

    if ((pLRate->mainDimension != pRRate->mainDimension) ||
	(pLRate->bgDimension != pRRate->bgDimension) ||
	(pLRate->numberOfElements != pRRate->numberOfElements))
    {
	Error("CopyRawDiffRateBound: inconsistent dimensions", PROGRAM_ERROR);
    }
    
    for (i = 0; i < pLRate->mainDimension; i++)
    {
	for (j = 0; j < pLRate->bgDimension; j++)
	{
	    pLRate->bound[i][j][0] = pRRate->bound[i][j][0];
	    pLRate->bound[i][j][1] = pRRate->bound[i][j][1];
	}
    }
}

void ClipRawDiffRate(RawDiffRate* pRate, const int newSize)
{
    int clip;
    int newNumberOfElements;
    RawDiffRate tempRate;
    int i;
    int j;
    int k;
    int counter;
    int newCounter;


    printf("Clipping RawDiffRate.  From %8i...  ", pRate->numberOfElements);

    clip = pRate->mainDimension - newSize;

    for (newNumberOfElements = 0, i = clip; i < pRate->mainDimension; i++)
    {
	for (j = 0; j < pRate->bgDimension; j++)
	{
	    if ((pRate->bound[i][j][0] != -1) && 
		(pRate->bound[i][j][1] >= clip))    
		/* not empty AND within range */
	    {
		int newLowerBound;
		newLowerBound = IMax(pRate->bound[i][j][0], clip);
		/* new lower bound cannot be smaller than the lower clip */
		newNumberOfElements += pRate->bound[i][j][1] -
		    newLowerBound + 1;
	    }
	}
    }
    /* determine new number of elements */
    
    NewRawDiffRate(&tempRate, newSize, pRate->bgDimension, 
        newNumberOfElements);

    for (i = 0; i < newSize; i++)
    {
	int iOld = i + clip;    /* (old) index for old rate */
	for (j = 0; j < pRate->bgDimension; j++)
	{
#ifdef DEBUG
	    CheckIndex(0, pRate->mainDimension, iOld, "ClipRawDiffRate");
#endif
	    if ((pRate->bound[iOld][j][0] != -1) && 
		(pRate->bound[iOld][j][1] >= clip))
	    {
		tempRate.bound[i][j][0] = 
		    IMax(pRate->bound[iOld][j][0], clip) - clip;
		tempRate.bound[i][j][1] = 
		    pRate->bound[iOld][j][1] - clip;
	    }
	    else
	    {
		tempRate.bound[i][j][0] = -1;
		tempRate.bound[i][j][1] = -1;
	    }
	}
    }
    /* set the bounds */

    /* fast-forward the index */
    for (counter = 0, newCounter = 0, i = 0; i < pRate->mainDimension; i++)
    {
	int iNew = i - clip;    /* (new) index for new rate */
	for (j = 0; j < pRate->bgDimension; j++)
	{
	    if (pRate->bound[i][j][0] != -1)
	    {
		for (k = pRate->bound[i][j][0]; 
		     k <= pRate->bound[i][j][1]; k++, counter++)
		{
		    int kNew = k - clip;    /* (new) k index */

#ifdef DEBUG
		    CheckIndex(0, pRate->numberOfElements, counter,
			       "ClipRawDiffRate");
#endif
		    /* i must be within clips */
		    if (i >= clip)
		    {
			/* its own bound must not be empty && 
			   k must be within range: 
			   comparison w/ upper bound is not necessary 
			   because it is automatically satisfied */
#ifdef DEBUG
			CheckIndex(0, tempRate.mainDimension, iNew,
				   "ClipRawDiffRate");
#endif
			if ((tempRate.bound[iNew][j][0] != -1) &&
			    (kNew >= tempRate.bound[iNew][j][0]))
			{
#ifdef DEBUG
			    CheckIndex(0, tempRate.numberOfElements,
				       newCounter, "ClipRawDiffRate");
#endif
			    tempRate.diffRate[newCounter] = 
			        pRate->diffRate[counter];
			    newCounter++;
			}
		    }
		}
	    }
	}
    }
    if ((counter != pRate->numberOfElements) || 
	(newCounter != tempRate.numberOfElements))
    {
	Error("ClipRawDiffRate: counting does not match", PROGRAM_ERROR);
    }
    /* redundant checking, but comes cheap so why not */

    CopyRawDiffRate(pRate, &tempRate);
    /* this is where real copy takes place */
    DeleteRawDiffRate(&tempRate);

    printf("to %8i.\n", pRate->numberOfElements);
}

void EnlargeRawDiffRate(RawDiffRate* pRate, const int newSize)
{
    int i;
    int j;
    int numAddedCells;
    RawDiffRate tempRate;

    numAddedCells = newSize - pRate->mainDimension;

    NewRawDiffRate(&tempRate, newSize, pRate->bgDimension,
        pRate->numberOfElements);

    for (i = 0; i < newSize; i++)
    {
	int iOld = i - numAddedCells;
	for (j = 0; j < pRate->bgDimension; j++)
	{
	    if (iOld >= 0)    /* area where the original tensor is */
	    {
		if (pRate->bound[iOld][j][0] != -1)
		{
#ifdef DEBUG
		    CheckIndex(0, pRate->mainDimension, iOld, 
			       "EnlargeRawDiffRate");
#endif
		    tempRate.bound[i][j][0] = pRate->bound[iOld][j][0] +
		        numAddedCells;
		    tempRate.bound[i][j][1] = pRate->bound[iOld][j][1] +
		        numAddedCells;
		}
		else
		{
		    tempRate.bound[i][j][0] = -1;
		    tempRate.bound[i][j][1] = -1;
		}
	    }
	    else    /* outside the original tensor */
	    {
		tempRate.bound[i][j][0] = -1;
		tempRate.bound[i][j][1] = -1;
	    }
	}
    }

    for (i = 0; i < pRate->numberOfElements; i++)
    {
	tempRate.diffRate[i] = pRate->diffRate[i];
	/* the rate itself does not change */
    }

    CopyRawDiffRate(pRate, &tempRate);
    DeleteRawDiffRate(&tempRate);
}

void ModifyRawDiffRate(RawDiffRate* pRate, const int newSize)
{
    if (pRate->mainDimension > newSize)
    {
	ClipRawDiffRate(pRate, newSize);
    }
    else if (pRate->mainDimension < newSize)
    {
	EnlargeRawDiffRate(pRate, newSize);
    }
    else
    {
	/* do nothing */
    }
}

void ReadRawDiffRate(RawDiffRate* pRate, const char* filename)
{
    FILE* file;

    file = SafeFOpen(filename, "r");
    binfread((pRate->bound)[0][0], sizeof(int), (pRate->mainDimension)*
        (pRate->bgDimension)*2, file);
    binfread(pRate->diffRate, sizeof(double), pRate->numberOfElements, file);
    fclose(file);

    CheckRawDiffRate(pRate);
}



void NewTotalRate(TotalRate* pRate, const int num_main_bins)
{
    pRate->mainDimension = num_main_bins;
    pRate->totalRate = New_dVector(num_main_bins);

    InitializeTotalRate(pRate);
}

void DeleteTotalRate(TotalRate* pRate)
{
    Delete_dVector(pRate->totalRate);
}

void InitializeTotalRate(TotalRate* pRate)
{
    int i;

    for (i = 0; i < pRate->mainDimension; i++)
    {
	pRate->totalRate[i] = 0.;
    }
}

void CopyTotalRate(TotalRate* pLRate, const TotalRate* pRRate)
{
    int i;

    pLRate->mainDimension = pRRate->mainDimension;

    Delete_dVector(pLRate->totalRate);
    pLRate->totalRate = New_dVector(pRRate->mainDimension);

    for (i = 0; i < pRRate->mainDimension; i++)
    {
	pLRate->totalRate[i] = pRRate->totalRate[i];
    }
}

void ClipTotalRate(TotalRate* pRate, const int newSize)
{
    int i;
    int clip;
    TotalRate tempRate;


    printf("Clipping TotalRate.  From %5i...  ", pRate->mainDimension);
    /* let user know */

    clip = pRate->mainDimension - newSize;
    NewTotalRate(&tempRate, newSize);

    for (i = 0; i < newSize; i++)
    {
	int iOld = i + clip;
#ifdef DEBUG
	CheckIndex(0, pRate->mainDimension, iOld, "ClipTotalRate");
#endif
	tempRate.totalRate[i] = pRate->totalRate[iOld];
    }

    CopyTotalRate(pRate, &tempRate);
    DeleteTotalRate(&tempRate);
    
    printf("to %5i.\n", pRate->mainDimension);
}
    
void EnlargeTotalRate(TotalRate* pRate, const int newSize)
{
    int i;
    int numAddedCells;
    TotalRate tempRate;

    numAddedCells = newSize - pRate->mainDimension;

    NewTotalRate(&tempRate, newSize);

    for (i = 0; i < newSize; i++)
    {
	int iOld = i - numAddedCells;
	if (iOld >= 0)
	{
#ifdef DEBUG
	    CheckIndex(0, pRate->mainDimension, iOld, "EnlargeTotalRate");
#endif
	    tempRate.totalRate[i] = pRate->totalRate[iOld];
	}
	else
	{
	    tempRate.totalRate[i] = 0.;
	}
    }

    CopyTotalRate(pRate, &tempRate);
    DeleteTotalRate(&tempRate);
}

void ModifyTotalRate(TotalRate* pRate, const int newSize)
{
    if (pRate->mainDimension > newSize)
    {
	ClipTotalRate(pRate, newSize);
    }
    else if (pRate->mainDimension < newSize)
    {
	EnlargeTotalRate(pRate, newSize);
    }
    else
    {
	/* do nothing */
    }
}

void ReadTotalRate(TotalRate* pRate, const char* filename)
{
    FILE* file;

    file = SafeFOpen(filename, "r");
    binfread(pRate->totalRate, sizeof(double), (pRate->mainDimension), file);
    fclose(file);
}



void NewDiffRate(DiffRate* pRate, const int num_main_bins)
{
    pRate->mainDimension = num_main_bins;
    pRate->bound = New_iMatrix(num_main_bins, 2);
    pRate->diffRate = New_dMatrix(num_main_bins, num_main_bins);

    InitializeDiffRate(pRate);
}


void DeleteDiffRate(DiffRate* pRate)
{
    Delete_dMatrix(pRate->diffRate);
    Delete_iMatrix(pRate->bound);
}


void InitializeDiffRate(DiffRate* pRate)
{
    int i;
    int j;

    for (i = 0; i < pRate->mainDimension; i++)
    {
	pRate->bound[i][0] = pRate->mainDimension - 1;
	pRate->bound[i][1] = 0;
	/* note that the invalidation here is different from RawDiffRate */
	for (j = 0; j < pRate->mainDimension; j++)
	{
	    pRate->diffRate[i][j] = 0.;
	}
    }
}

void CopyDiffRate(DiffRate* pLRate, const DiffRate* pRRate)
{
    int i;
    int j;

    pLRate->mainDimension = pRRate->mainDimension;

    Delete_dMatrix(pLRate->diffRate);
    Delete_iMatrix(pLRate->bound);

    pLRate->bound = New_iMatrix(pRRate->mainDimension, 2);
    pLRate->diffRate = New_dMatrix(pRRate->mainDimension, 
        pRRate->mainDimension);

   for (i = 0; i < pRRate->mainDimension; i++)
    {
	pLRate->bound[i][0] = pRRate->bound[i][0];
	pLRate->bound[i][1] = pRRate->bound[i][1];
	for (j = 0; j < pRRate->mainDimension; j++)
	{
	    pLRate->diffRate[i][j] = pRRate->diffRate[i][j];
	}
    }
}

void CopyDiffRateBound(DiffRate* pLRate, const DiffRate* pRRate)
{
    int i;

    if (pLRate->mainDimension != pRRate->mainDimension)
    {
	Error("CopyDiffRateBound: inconsistent dimensions", PROGRAM_ERROR);
    }

    for (i = 0; i < pLRate->mainDimension; i++)
    {
	pLRate->bound[i][0] = pRRate->bound[i][0];
	pLRate->bound[i][1] = pRRate->bound[i][1];
    }
}

void SetStandardDiffRateBound(DiffRate* pRate)
{
    int i;
    
    InvalidateDiffRateBound(pRate);

    for (i = 0; i < pRate->mainDimension; i++)
    {
	if ((pRate->bound[i][0] != pRate->mainDimension - 1) ||
	    (pRate->bound[i][1] != 0))
	{
	    pRate->bound[i][0] = 0;
	    pRate->bound[i][1] = i;
	}
    }
}

void ClipDiffRate(DiffRate* pRate, const int newSize)
{
    int clip;
    DiffRate tempRate;
    int i;
    int j;

    printf("Clipping DiffRate.  From %5i...  ", pRate->mainDimension);

    clip = pRate->mainDimension - newSize;

    NewDiffRate(&tempRate, newSize);

    for (i = 0; i < newSize; i++)
    {
	int iOld = i + clip;    /* (old) index for old rate */

#ifdef DEBUG
	CheckIndex(0, pRate->mainDimension, iOld, "ClipDiffRate");
#endif
	if ((pRate->bound[iOld][0] != pRate->mainDimension - 1) ||
	    (pRate->bound[iOld][1] != 0))
	/* bound is valid */
	{
	    if (pRate->bound[iOld][1] >= clip)
	    {
		tempRate.bound[i][0] = IMax(pRate->bound[iOld][0], clip) -
		    clip;
		tempRate.bound[i][1] = pRate->bound[iOld][1] - clip;
	    }
	    else
	    {
		tempRate.bound[i][0] = newSize - 1;
		tempRate.bound[i][1] = 0;
	    }
	}
	else    /* bound is invalid to begin with */
	{
	    tempRate.bound[i][0] = newSize - 1;
	    tempRate.bound[i][1] = 0;
	}

	for (j = 0; j < newSize; j++)
	{
	    int jOld = j + clip;    /* (old) j-index */

#ifdef DEBUG
	    CheckIndex(0, pRate->mainDimension, jOld, "ClipDiffRate");
#endif
	    tempRate.diffRate[i][j] = pRate->diffRate[iOld][jOld];
	}
    }

    CopyDiffRate(pRate, &tempRate);
    DeleteDiffRate(&tempRate);

    printf("to %5i.\n", pRate->mainDimension);
}

void EnlargeDiffRate(DiffRate* pRate, const int newSize)
{
    int i;
    int j;
    int numAddedCells;
    DiffRate tempRate;

    numAddedCells = newSize - pRate->mainDimension;

    NewDiffRate(&tempRate, newSize);

    for (i = 0; i < newSize; i++)
    {
	int iOld = i - numAddedCells;
	if (iOld >= 0)
	{
#ifdef DEBUG
	    CheckIndex(0, pRate->mainDimension, iOld, "EnlargeDiffRate");
#endif
	    if ((pRate->bound[iOld][0] != pRate->mainDimension - 1) ||
		(pRate->bound[iOld][1] != 0))
	    /* old bound is valid */
	    {
		tempRate.bound[i][0] = pRate->bound[iOld][0] + numAddedCells;
		tempRate.bound[i][1] = pRate->bound[iOld][1] + numAddedCells;
	    }
	    else    /* old bound is invalild */
	    {
		tempRate.bound[i][0] = newSize - 1;
		tempRate.bound[i][1] = 0;
	    }
	}
	else
	{
	    tempRate.bound[i][0] = newSize - 1;
	    tempRate.bound[i][1] = 0;
	}

	for (j = 0; j < newSize; j++)
	{
	    int jOld = j - numAddedCells;
	    if (iOld >= 0 && jOld >= 0)
	    {
#ifdef DEBUG
		CheckIndex(0, pRate->mainDimension, jOld, "EnlargeDiffRate");
#endif
		tempRate.diffRate[i][j] = pRate->diffRate[iOld][jOld];
	    }
	    else
	    {
		tempRate.diffRate[i][j] = 0.;
	    }
	}
    }

    CopyDiffRate(pRate, &tempRate);
    DeleteDiffRate(&tempRate);
}

void ModifyDiffRate(DiffRate* pRate, const int newSize)
{
    if (pRate->mainDimension > newSize)
    {
	ClipDiffRate(pRate, newSize);
    }
    else if (pRate->mainDimension < newSize)
    {
	EnlargeDiffRate(pRate, newSize);
    }
    else
    {
	/* do nothing */
    }
}

void ReadDiffRate(DiffRate* pRate, const char* filename)
{
    FILE* file;
    
    file = SafeFOpen(filename, "r");
    binfread((pRate->bound)[0], sizeof(int), (pRate->mainDimension)*2, file);
    binfread((pRate->diffRate)[0], sizeof(double), (pRate->mainDimension)*
        (pRate->mainDimension), file);
    fclose(file);

    InvalidateDiffRateBound(pRate);
    /* make sure invalid bounds are properly invalidated */
}

void InvalidateDiffRateBound(DiffRate* pRate)
{
    int i;

    for (i = 0; i < pRate->mainDimension; i++)
    {
	if ((pRate->bound[i][0] == -1) || ((pRate->bound[i][0] ==
	    pRate->mainDimension - 1) && (pRate->bound[i][0] == 0)))
	/* picks up invalid bound data however they are formatted */
	{
	    pRate->bound[i][0] = pRate->mainDimension - 1;
	    pRate->bound[i][1] = 0;
	}
    }
}


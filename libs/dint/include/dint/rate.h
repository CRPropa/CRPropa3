#ifndef _RATE_H
#define _RATE_H

#include "dint/vector.h"
#include "dint/binfread.h"

typedef struct
{
    int mainDimension;
    int bgDimension;
    //    int mainDimension2;    /* number of daughter bins */
    int numberOfElements;
    i3Tensor bound;
    dVector diffRate;
} RawDiffRate;
/* although RawDiffRate is really a 3-dimensional matrix, I "serialize"
    the matrix to tighten up the memory: instead of using the full 
    3-dimensional matrix with lots of zeros, I replace it by a serialized
    1-dimensional array with non-zero elements and with appropriate lower
    and upper bounds (bound) */

typedef struct
{
    int mainDimension;
    int bgDimension;
    dMatrix totalRate;
} RawTotalRate;

typedef struct
{
    int mainDimension;
    iMatrix bound;
    dMatrix diffRate;
} DiffRate;

typedef struct
{
    int mainDimension;
    dVector totalRate;
} TotalRate;
/* NOTE: all rates are initialized when assigned memory.  All bounds are
   properly invalidated. */


void NewRawTotalRate(RawTotalRate* pRate, const int num_main_bins,
		     const int num_bg_bins);
void DeleteRawTotalRate(RawTotalRate* pRate);
void InitializeRawTotalRate(RawTotalRate* pRate);
void CopyRawTotalRate(RawTotalRate* pLRate, const RawTotalRate* pRRate);
void ClipRawTotalRate(RawTotalRate* pRate, const int newSize);
void EnlargeRawTotalRate(RawTotalRate* pRate, const int newSize);
void ModifyRawTotalRate(RawTotalRate* pRate, const int newSize);
void ReadRawTotalRate(RawTotalRate* pRate, const char* filename);


void NewRawDiffRate(RawDiffRate* pRate, const int num_main_bins, 
		    const int num_bg_bins, const int num_elements);
void DeleteRawDiffRate(RawDiffRate* pRate);
void InitializeRawDiffRate(RawDiffRate* pRate);
void CheckRawDiffRate(RawDiffRate* pRate);
#ifdef DEBUG
double RawDiffRateElement(const RawDiffRate* pRate, const int i, const int j,
			  const int k);
/* although this is a very good and safe way of getting the element
   it is a very expensive call: it is recommended only for checking */
#endif
void CopyRawDiffRate(RawDiffRate* pLRate, const RawDiffRate* pRRate);
void CopyRawDiffRateBound(RawDiffRate* pLRate, const RawDiffRate* pRRate);
void ClipRawDiffRate(RawDiffRate* pRate, const int newSize);
void EnlargeRawDiffRate(RawDiffRate* pRate, const int newSize);
void ModifyRawDiffRate(RawDiffRate* pRate, const int newSize);
void ReadRawDiffRate(RawDiffRate* pRate, const char* filename);


void NewTotalRate(TotalRate* pRate, const int num_main_bins);
void DeleteTotalRate(TotalRate* pRate);
void InitializeTotalRate(TotalRate* pRate);
void CopyTotalRate(TotalRate* pLRate, const TotalRate* pRRate);
void ClipTotalRate(TotalRate* pRate, const int newSize);
void EnlargeTotalRate(TotalRate* pRate, const int newSize);
void ModifyTotalRate(TotalRate* pRate, const int newSize);
void ReadTotalRate(TotalRate* pRate, const char* filename);


void NewDiffRate(DiffRate* pRate, const int num_main_bins);
void DeleteDiffRate(DiffRate* pRate);
void InitializeDiffRate(DiffRate* pRate);
void CopyDiffRate(DiffRate* pLRate, const DiffRate* pRRate);
void ClipDiffRate(DiffRate* pRate, const int newSize);
void EnlargeDiffRate(DiffRate* pRate, const int newSize);
void ModifyDiffRate(DiffRate* pRate, const int newSize);
void ReadDiffRate(DiffRate* pRate, const char* filename);
void CopyDiffRateBound(DiffRate* pLRate, const DiffRate* pRRate);
void SetStandardDiffRateBound(DiffRate* pRate);
/* sets DiffRate bound to 0 <= k <= i */
void InvalidateDiffRateBound(DiffRate* pRate);

#endif
/* Copy*Rate: perform real deep copies; i.e. they reformat the dimensions
   according to the size of the rate being copied. 
   Copy*DiffRateBound: copy bounds without altering size.
*/

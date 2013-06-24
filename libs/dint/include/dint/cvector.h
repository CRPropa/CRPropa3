
#ifndef _CVECTOR_H
#define _CVECTOR_H

#include "dint/vector.h"
#include <stdio.h>

// these are variants of structs defined in vector.h; they have dimensional
//   info built-in, and they are initialized when assigned memory
typedef struct
{
    int dimension;
    dVector vector;
} dCVector;

typedef struct
{
    int dimension1;
    int dimension2;
    iMatrix matrix;
} iCMatrix;

void New_dCVector(dCVector* pVector, const int n);
void Delete_dCVector(dCVector* pVector);
void Initialize_dCVector(dCVector* pVector);

void New_iCMatrix(iCMatrix* pMatrix, const int n1, const int n2);
void Delete_iCMatrix(iCMatrix* pMatrix);
void Initialize_iCMatrix(iCMatrix* pMatrix);

#endif

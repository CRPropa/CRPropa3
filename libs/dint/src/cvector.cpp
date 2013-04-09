
#include "dint/cvector.h"

void New_dCVector(dCVector* pVector, const int n) {
  pVector->dimension = n;
  pVector->vector = New_dVector(n);
  
  Initialize_dCVector(pVector);
}

void Delete_dCVector(dCVector* pVector) {
  Delete_dVector(pVector->vector);
}

void Initialize_dCVector(dCVector* pVector) {
  for (int i = 0; i<pVector->dimension; i++) pVector->vector[i] = 0.;
}


void New_iCMatrix(iCMatrix* pMatrix, const int n1, const int n2) {
  pMatrix->dimension1 = n1;
  pMatrix->dimension2 = n2;
  pMatrix->matrix = New_iMatrix(n1, n2);
  
  Initialize_iCMatrix(pMatrix);
}

void Delete_iCMatrix(iCMatrix* pMatrix) {
  Delete_iMatrix(pMatrix->matrix);
}

void Initialize_iCMatrix(iCMatrix* pMatrix) {
  for (int i=0; i<pMatrix->dimension1; i++) {
    for (int j=0; j<pMatrix->dimension2; j++) pMatrix->matrix[i][j] = 0;
  }
}

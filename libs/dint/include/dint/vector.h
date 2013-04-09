#ifndef _VECTOR_H_
#define _VECTOR_H_

/* Vectors, matrices, and tensors defined here are very simple, dynamic
   objects.  They do not provide any bound checking nor has extra data
   structure, but have very little overhead because of that.  For structs
   that have dimensional data built-in, see "cvector.h".
*/

typedef double* dVector;
typedef int* iVector;
typedef double** dMatrix;
typedef int** iMatrix;
typedef double*** d3Tensor;
typedef int*** i3Tensor;


dVector New_dVector(const int n);
iVector New_iVector(const int n);
dMatrix New_dMatrix(const int n1, const int n2);
iMatrix New_iMatrix(const int n1, const int n2);
d3Tensor New_d3Tensor(const int n1, const int n2, const int n3);
i3Tensor New_i3Tensor(const int n1, const int n2, const int n3);
void Delete_dVector(dVector vector);
void Delete_iVector(iVector vector);
void Delete_dMatrix(dMatrix matrix);
void Delete_iMatrix(iMatrix matrix);
void Delete_d3Tensor(d3Tensor tensor);
void Delete_i3Tensor(i3Tensor tensor);

#endif

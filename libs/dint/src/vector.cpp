
#include <stdio.h>
#include <stdlib.h>

#include "dint/vector.h"
#include "dint/error.h"


dVector New_dVector(const int n)
{
    dVector vector;

    vector = (double*)malloc((size_t)(n*sizeof(double)));
    if (vector == NULL)
    {
	Error("New_dVector: allocation failure", ARRAY_ERROR);
    }

    return vector;
}

iVector New_iVector(const int n)
{
    iVector vector;

    vector = (int*)malloc((size_t)(n*sizeof(int)));
    if (vector == NULL)
    {
	Error("New_iVector: allocation failure", ARRAY_ERROR);
    }

    return vector;
}


dMatrix New_dMatrix(const int n1, const int n2)
{
    int i;
    dMatrix matrix;

    matrix = (double**)malloc((size_t)(n1*sizeof(double*)));
    if (matrix == NULL)
    {
	Error("New_dMatrix: allocation failure level 1", ARRAY_ERROR);
    }

    matrix[0] = (double*)malloc((size_t)(n1*n2*sizeof(double)));
    if (matrix[0] == NULL)
    {
	Error("New_dMatrix: allocation failure level 2", ARRAY_ERROR);
    }
	
    for (i = 1; i < n1; i++)
    {
        matrix[i] = matrix[i-1] + n2;
    }

    return matrix;
}


iMatrix New_iMatrix(const int n1, const int n2)
{
    int i;
    iMatrix matrix;

    matrix = (int**)malloc((size_t)(n1*sizeof(int*)));
    if (matrix == NULL)
    {
	Error("New_iMatrix: allocation failure level 1", ARRAY_ERROR);
    }

    matrix[0] = (int*)malloc((size_t)(n1*n2*sizeof(int)));
    if (matrix[0] == NULL)
    {
	Error("New_iMatrix: allocation failure level 2", ARRAY_ERROR);
    }

    for (i = 1; i < n1; i++)
    {
        matrix[i] = matrix[i-1] + n2;
    }

    return matrix;
}


d3Tensor New_d3Tensor(const int n1, const int n2, const int n3)
{
    int i;
    int j;
    d3Tensor tensor;

    tensor = (double***)malloc((size_t)(n1*sizeof(double**)));
    if (tensor == NULL)
    {
	Error("d3Tensor: allocation failure level 1", ARRAY_ERROR);
    }

    tensor[0] = (double**)malloc((size_t)(n1*n2*sizeof(double*)));
    if (tensor[0] == NULL)
    {
	Error("d3Tensor: allocation failure level 2", ARRAY_ERROR);
    }

    tensor[0][0] = (double*)malloc((size_t)(n1*n2*n3*sizeof(double)));
    if (tensor[0][0] == NULL)
    {
	Error("d3Tensor: allocation failure level 3", ARRAY_ERROR);
    }

    for (j = 1; j < n2; j++)
    {
        tensor[0][j] = tensor[0][j-1] + n3;
    }
    for (i = 1; i < n1; i++)
    {
        tensor[i] = tensor[i-1] + n2;
	tensor[i][0] = tensor[i-1][0] + n2*n3;
	for (j = 1; j < n2; j++)
	{
	    tensor[i][j] = tensor[i][j-1] + n3;
	}
    }

    return tensor;
}


i3Tensor New_i3Tensor(const int n1, const int n2, const int n3)
{
    int i;
    int j;
    i3Tensor tensor;

    tensor = (int***)malloc((size_t)(n1*sizeof(int**)));
    if (tensor == NULL)
    {
	Error("New_i3Tensor: allocation failure level 1", ARRAY_ERROR);
    }

    tensor[0] = (int**)malloc((size_t)(n1*n2*sizeof(int*)));
    if (tensor[0] == NULL)
    {
	Error("New_i3Tensor: allocation failure level 2", ARRAY_ERROR);
    }

    tensor[0][0] = (int*)malloc((size_t)(n1*n2*n3*sizeof(int)));
    if (tensor[0][0] == NULL)
    {
	Error("New_i3Tensor: allocation failure level 3", ARRAY_ERROR);
    }


    for (j = 1; j < n2; j++)
    {
        tensor[0][j] = tensor[0][j-1] + n3;
    }
    for (i = 1; i < n1; i++)
    {
        tensor[i] = tensor[i-1] + n2;
	tensor[i][0] = tensor[i-1][0] + n2*n3;
	for (j = 1; j < n2; j++)
	{
	    tensor[i][j] = tensor[i][j-1] + n3;
	}
    }

    return tensor;
}


void Delete_dVector(dVector vector)
{
    free(vector);
    vector = NULL;
}


void Delete_iVector(iVector vector)
{
    free(vector);
    vector = NULL;
}


void Delete_dMatrix(dMatrix matrix)
{
    free(matrix[0]);
    matrix[0] = NULL;
    free(matrix);
    matrix = NULL;
}


void Delete_iMatrix(iMatrix matrix)
{
    free(matrix[0]);
    matrix[0] = NULL;
    free(matrix);
    matrix = NULL;
}


void Delete_d3Tensor(d3Tensor tensor)
{
    free(tensor[0][0]);
    tensor[0][0] = NULL;
    free(tensor[0]);
    tensor[0] = NULL;
    free(tensor);
    tensor = NULL;
}


void Delete_i3Tensor(i3Tensor tensor)
{
    free(tensor[0][0]);
    tensor[0][0] = NULL;
    free(tensor[0]);
    tensor[0] = NULL;
    free(tensor);
    tensor = NULL;
}

#ifndef __UTILS_H__
#define __UTILS_H__

#include <iostream>
#include <cblas.h>
#include <cstring>
#include <math.h>

enum MatOrder
{
    colMajor = 0,
    rowMajor = 1
};

template<typename T>
inline void printVec(T* vec, size_t size)
{
    std::cout << "[ ";
    for(int i=0; i<size-1; i++)
        std::cout << vec[i] << " ";
    if(size-1 >= 0)
        std::cout << vec[size-1];
    std::cout << " ]\n";
}

template<typename T>
inline void printMat(T* mat, size_t rowSize, size_t colSize, MatOrder order)
{

    if(order == MatOrder::rowMajor)
    {
        for(size_t i=0; i<rowSize; i++)
        {
            for(size_t j=0; j<colSize; j++)
            {
                std::cout << mat[i*colSize + j] << " ";
            }
            std::cout << "\n";
        }
    }
    else if(order == MatOrder::colMajor)
    {
        for(size_t i=0; i<rowSize; i++)
        {
            for(size_t j=0; j<colSize; j++)
            {
                std::cout << mat[j*colSize + i] << " ";
            }
            std::cout << "\n";
        }
    }
}


/**
 * Create Identity matrix by the specified size
 * 
 * Parameters :
 * size
 * 
 * Returns:
 * 
 * resulting identity matrix
 */
double* createIdentity(int size)
{
    double* ans = (double*)calloc(size*size, sizeof(double));

    for(int i=0; i<size; i++)
        ans[i*size + i] = 1.0;

    return ans;
}

/**
 * Returns a zero vector
 * 
 * Parameters:
 * 
 * size -> size
 * 
 * Returns:
 * 
 * double* vector
 */
inline double* zeros(size_t size)
{
    return (double*)calloc(size,sizeof(double));
}

/**
 * Returns a zero matrix
 * 
 * Parameters:
 * rowSize    -> row size
 * columnSize -> column size
 * 
 * Returns:
 * 
 * double* matrix
 */
inline double* zeros(size_t rowSize, size_t columnSize)
{
    return (double*)calloc(rowSize*columnSize, sizeof(double));
}


inline double* ones(size_t size)
{  
    double* ans = (double*)malloc(sizeof(double)*size);

    for(size_t i=0; i<size; i++)
        ans[i] = 1.0;

    return ans;
}

inline double* ones(size_t rowSize, size_t columnSize)
{
    double* ans = (double*)malloc(sizeof(double)*rowSize*columnSize);

    for(size_t i=0; i<rowSize*columnSize; i++)
        ans[i] = 1.0;

    return ans;
}

/**
 * 
 * Gets lower triangular part of matrix
 * 
 * Parameters:
 * 
 * A      -> matrix
 * size   -> matrix size
 * offset -> determines distance from diagonal
 * 
 * Returns:
 * 
 * double *ans -> result matrix
 * 
 */
inline double* tril(double *A, int size, int offset)
{
    double* ans = zeros(size, size);
    for(int i=0; i<size; i++)
    {
        for(int j=i-offset; j>=0; j--)
            ans[i*size + j] = A[i*size + j];
    }

    return ans;
}

/**
 * 
 * Gets upper triangular part of matrix
 * 
 * Parameters:
 * 
 * A      -> matrix
 * size   -> matrix size
 * offset -> determines distance from diagonal
 * 
 * Returns:
 * 
 * double *ans -> result matrix
 * 
 */
inline double* triu(double *A, size_t size, int offset)
{
    double* ans = zeros(size, size);
    for(size_t i=0; i<size; i++)
    {
        for(size_t j=i+offset; j<size; j++)
            ans[i*size + j] = A[i*size + j];
    }

    return ans;
}

/**
 * 
 * Gets diagonal part of matrix
 * 
 * Parameters:
 * 
 * A      -> matrix
 * size   -> matrix size
 * 
 * Returns:
 * 
 * double *ans -> result matrix
 * 
 */
inline double *diag(double *A, size_t size)
{
    double* ans = zeros(size, size);

    for(size_t i=0; i<size; i++)
        ans[i*size + i] = A[i*size + i];

    return ans;
}

/**
 * 
 * Change Order of matrix from Column major to Row major 
 * 
 */
inline void changeCol2Row(double * A, int rowSize, int colSize)
{
    double *changedOrderMatrix = zeros(rowSize, colSize);

    for(int i=0; i<rowSize; i++)
    {
        for(int j=0; j<colSize; j++)
            changedOrderMatrix[i*colSize + j] = A[j*rowSize + i];
    }

    std::memcpy(A, changedOrderMatrix, sizeof(double)*rowSize*colSize);
    free(changedOrderMatrix);
}

#endif
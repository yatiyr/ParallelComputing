#ifndef __SEQUENTIAL_FUNCTIONS_H__
#define __SEQUENTIAL_FUNCTIONS_H__

#include <Utils.h>

/**
 * 
 * Adds matrices
 * 
 * Parameters:
 * 
 * matrix1
 * matrix2
 * rowSize
 * columnSize
 * 
 * Returns
 * 
 * writes result to matrix2
 * 
 */
inline void addMatrices_seq(double* matrix1, double* matrix2, int rowSize, int columnSize)
{
    for(int i=0; i<rowSize; i++)
    {
        for(int j=0; j<columnSize; j++)
        {
            matrix2[i*columnSize + j] += matrix1[i*columnSize + j];
        }
    }
}

/**
 * This function performs forward sweep operation on a lower triangular matrix
 * Matrices are assumed in row major order
 * 
 * It is assumed that memory is allocated for vector x
 * 
 * Parameters:
 * 
 * m     -> lower triangular matrix
 * mSize -> size of the matrix
 * b     -> right hand side vector of system
 * x     -> result vector
 */
inline void forward_sweep_seq(double *A, int mSize, double *b, double *x)
{
    // set x values to 0
    std::memset(x, 0, sizeof(double)*mSize);
    
    for(int i=0; i<mSize; i++)
    {
        for(int j=i-1; j>=0; j--)
            x[i] += A[i * mSize + j]*x[j];

        x[i] = (b[i] - x[i])/A[i* mSize + i];
    }
}

/**
 * Solves a linear system Ax = b with gauss-seidel method sequentially
 * 
 * It is assumed that matrices are in row major order
 * 
 * Used CBLAS Functions :
 * 
 * dgemv -> for matrix vector mult
 * daxpy -> for taking difference between vectors
 * dnrm2 -> get norm of the vector
 * 
 * Parameters:
 * 
 * A             -> Matrix (diagonally dominant and spectral radius is < 1 for ensuring convergence)
 * b             -> right hand side vector of system
 * x0            -> initial guess vector
 * ans           -> answer vector
 * size          -> size of matrix
 * maxIterations -> maximum number of iterations
 * tolerance     -> error tolerance
 * Returns:
 * 
 * double *x0         -> resulting vector
 * int iterationCount -> how many iterations passed
 * 
 * 
 */
inline int gauss_seidel_seq(double *A, double *b, double *x0, size_t size, size_t maxIterations, double tolerance)
{

    // extract lower triangular part
    // diagonal is included
    double* L = tril(A, size, 0);

    // extract upper triangular part
    // diagonal is excluded
    double* U = triu(A, size, 1);

    // allocate memory for x and fill it with 
    // initial guess
    double* x = (double*) malloc(sizeof(double)*size);
    std::memcpy(x, x0, sizeof(double)*size);

    // allocate memory for x in next iterations
    // to use inside for loop
    double* x_new = (double*) malloc(sizeof(double)*size);

    // allocate memory for b_new and fill it with
    // our right hand side vector b
    double* b_new = (double*) malloc(sizeof(double)*size);

    // allocate memory for difference vector
    double* diff = (double*) malloc(sizeof(double)*size);


    // this keeps track of how many iterations have passed
    int iterations = 0;

    for(size_t i=0; i<maxIterations; i++)
    {
        // get b_new equal to b again
        std::memcpy(b_new, b, sizeof(double)*size);        

        // This cblas call calculates (b - U*x) and writes the result to b_new
        cblas_dgemv(CBLAS_ORDER::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans, size, size, -1.0, U, size, x, 1, 1.0, b_new, 1);

        // get new x vector
        forward_sweep_seq(L, size, b_new, x_new);

        // increase iteration count
        iterations++;
        
        // set elements of x_new to diff vector
        std::memcpy(diff, x_new, sizeof(double)*size);

        // get difference vector between x and x_new
        // and write it to diff vector
        cblas_daxpy(size, -1, x, 1, diff, 1);

        // if norm of diff smaller than tolerance
        // update x with x_new and get out of the loop
        if(cblas_dnrm2(size, diff, 1) < tolerance)
        {
            std::memcpy(x, x_new, sizeof(double)*size);
            break;
        }

        // update x with x_new
        std::memcpy(x, x_new, sizeof(double)*size);

    }

    // copy result to x0 array
    std::memcpy(x0, x, sizeof(double)*size);

    // deallocate tmp arrays
    free(L);
    free(U);
    free(x);
    free(x_new);
    free(b_new);
    free(diff);

    return iterations;
}

/**
 * Solves a linear system Ax = b with jacobi method sequentially
 * 
 * It is assumed that matrices are in row major order
 * 
 * Used CBLAS Functions :
 * 
 * dgemv -> for matrix vector mult
 * daxpy -> for taking difference between vectors
 * dnrm2 -> get norm of the vector
 * 
 * Parameters:
 * 
 * A             -> Matrix (diagonally dominant and spectral radius is < 1 for ensuring convergence)
 * b             -> right hand side vector of system
 * x0            -> initial guess vector
 * ans           -> answer vector
 * size          -> size of matrix
 * maxIterations -> maximum number of iterations
 * tolerance     -> error tolerance
 * Returns:
 * 
 * double *x0         -> resulting vector
 * int iterationCount -> how many iterations passed
 * 
 * 
 */
inline int jacobi_seq(double *A, double *b, double *x0, size_t size, size_t maxIterations, double tolerance)
{

    // extract lower triangular part
    // diagonal is excluded
    double* L = tril(A, size, 1);

    // extract upper triangular part
    // diagonal is excluded
    double* U = triu(A, size, 1);

    // extract diagonal part of the
    // matrix
    double* D = diag(A, size);

    // add Lower and Upper parts together
    // and store result in U
    addMatrices_seq(L,U,size,size);

    // allocate memory for x and fill it with 
    // initial guess
    double* x = (double*) malloc(sizeof(double)*size);
    std::memcpy(x, x0, sizeof(double)*size);

    // allocate memory for x in next iterations
    // to use inside for loop
    double* x_new = (double*) malloc(sizeof(double)*size);

    // allocate memory for b_new and fill it with
    // our right hand side vector b
    double* b_new = (double*) malloc(sizeof(double)*size);

    // allocate memory for difference vector
    double* diff = (double*) malloc(sizeof(double)*size);


    // this keeps track of how many iterations have passed
    int iterations = 0;

    for(size_t i=0; i<maxIterations; i++)
    {
        // get b_new equal to b again
        std::memcpy(b_new, b, sizeof(double)*size);        

        // This cblas call calculates (b - U*x) and writes the result to b_new
        cblas_dgemv(CBLAS_ORDER::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans, size, size, -1.0, U, size, x, 1, 1.0, b_new, 1);

        // get new x vector
        forward_sweep_seq(D, size, b_new, x_new);

        // increase iteration count
        iterations++;
        
        // set elements of x_new to diff vector
        std::memcpy(diff, x_new, sizeof(double)*size);

        // get difference vector between x and x_new
        // and write it to diff vector
        cblas_daxpy(size, -1, x, 1, diff, 1);

        // if norm of diff smaller than tolerance
        // update x with x_new and get out of the loop
        if(cblas_dnrm2(size, diff, 1) < tolerance)
        {
            std::memcpy(x, x_new, sizeof(double)*size);
            break;
        }

        // update x with x_new
        std::memcpy(x, x_new, sizeof(double)*size);

    }

    // copy result to x0 array
    std::memcpy(x0, x, sizeof(double)*size);

    // deallocate tmp arrays
    free(L);
    free(U);
    free(D);
    free(x);
    free(x_new);
    free(b_new);
    free(diff);

    return iterations;
}

/**
 * Solves a linear system Ax = b with  chebyshev iteration sequentially
 * 
 * It is assumed that matrices are in row major order
 * 
 * Parameters:
 * 
 * A             -> Matrix (diagonally dominant and spectral radius is < 1 for ensuring convergence)
 * b             -> right hand side vector of system
 * x0            -> initial guess vector
 * ans           -> answer vector
 * size          -> size of matrix
 * maxIterations -> maximum number of iterations
 * tolerance     -> error tolerance
 * lMax          -> largest estimated eigenvalue of matrix A
 * lMin          -> smallest estimated eigenvalue of Matrix A
 * Returns:
 * 
 * double *x0         -> resulting vector
 * int iterationCount -> how many iterations passed
 * 
 * 
 */
inline int chebyshev_seq(double *A, double *b, double *x0, size_t size, int maxIterations, double tolerance, double lMax, double lMin)
{


    double d = (lMax + lMin) / 2;
    double c = (lMax - lMin) / 2;

    double alpha = 0;
    double beta  = 0;

    // preconditioner is an identity matrix
    // of size(A)
    double* preCond = createIdentity(size);

    // allocate memory for x and fill it with 
    // initial guess
    double* x = (double*) malloc(sizeof(double)*size);
    std::memcpy(x, x0, sizeof(double)*size);

    // allocate memory for x in next iterations
    // to use inside for loop
    double* x_new = (double*) malloc(sizeof(double)*size);

    // allocate memory for difference vector
    double* diff = (double*) malloc(sizeof(double)*size);

    // allocate memory for z vector
    double* z = (double*) malloc(sizeof(double)*size);

    // allocate memory for p vector
    double* p = (double*) malloc(sizeof(double)*size);

    // allocate memory for initial guess of precondition
    // system
    double* x_pre = zeros(size);
    double* x_preNew = (double*) malloc(sizeof(double)*size);

    // allocate value for r
    double* r = (double*) malloc(sizeof(double)*size);
    std::memcpy(r,b,sizeof(double)*size);

    // This cblas call calculates (r = b - A*x) and writes the result to r
    cblas_dgemv(CBLAS_ORDER::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans, size, size, -1.0, A, size, x, 1, 1.0, r, 1);    

    // this keeps track of how many iterations have passed
    int iterations = 0;

    for(int i=0; i<maxIterations; i++)
    {
        iterations++;

        // fill x_pre with zeros
        std::memcpy(x_preNew,x_pre,sizeof(double)*size);
        gauss_seidel_seq(preCond, r, x_preNew, size, 100, 1e-15);

        if(i == 0)
        {
            std::memcpy(p, x_preNew, sizeof(double)*size);
            alpha = 1/d;
        }
        else if(i == 1)
        {
            beta  = 0.5 * std::pow(c * alpha, 2);
            alpha =  1 / (d - beta / alpha);

            // This cblas call calculates (p = z + beta*p) and writes the result to p
            cblas_daxpy(size, beta, p, 1, p, 1);            
        }
        else
        {
            beta  = std::pow(c*alpha/2,2);
            alpha = 1/(d - beta/alpha);

            // These cblas calls calculate (p = z + beta*p) and writes the result to p
            cblas_dscal(size, beta, p , 1);
            cblas_daxpy(size, 1, x_preNew, 1, p, 1);
        }

        // x = x + alpha * p
        cblas_daxpy(size, alpha, p, 1, x, 1);

        // make r equal to b again
        std::memcpy(r, b, sizeof(double)*size);

        // r = b - A*x
        cblas_dgemv(CBLAS_ORDER::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans, size, size, -1.0, A, size, x, 1, 1, r, 1);
        double norm_r = cblas_dnrm2(size, r, 1);
        if(norm_r < tolerance)
        {
            break;
        }

    }

    std::memcpy(x0, x, sizeof(double)*size);


    free(preCond);
    free(x);
    free(x_new);
    free(diff);
    free(z);
    free(p);
    free(x_pre);
    free(x_preNew);
    free(r);

    return iterations;
}


#endif
#ifndef __PARALLEL_FUNCTIONS_H__
#define __PARALLEL_FUNCTIONS_H__

#include <mpi.h>
#include <cblas.h>
#include <SequentialFunctions.h>

#define MANAGER_RANK 0

#define MANAGER_TAG 1
#define WORKER_TAG 2

// Parallel_cblas_dgemv
// OpenMPI implementation of matrix vector multiplication subrouite of blas
// It takes ordinary dgemv parameters but we also pass offset and chunksize
// for partitioning and we gather result with allgather
// We do row partitioning for the matrix
inline void parallel_cblas_dgemv(CBLAS_ORDER order, CBLAS_TRANSPOSE trans, size_t M, size_t N, 
                                 double Alpha, double* A, size_t lda, double* X,
                                 int incX, double Beta, double* Y, int incY, int offset, int chunksize, int colSize)
{
    // Calculate offset of matrix
    int matrixOffset    = offset * colSize;

    // Get pointers of partitions
    double* matrixPart  = &A[matrixOffset];
    double* xVectorPart = &X[offset];
    double* yVectorPart = &Y[offset];

    // Calculates matrix vector operation and writes the result to corresponding part
    // of Y vector
    cblas_dgemv(order, trans, chunksize, colSize, Alpha, matrixPart, colSize, xVectorPart,
                incX, Beta, yVectorPart, incY);

    // We gather all parts of Y which is calculated in every process and write it into Y itself
    MPI_Allgather(yVectorPart, chunksize, MPI_DOUBLE, Y, chunksize, MPI_DOUBLE, MPI_COMM_WORLD);

}


// Parallel cblas dscal
// performs dscal operation in parallel
inline void parallel_cblas_dscal(int Size, double Alpha, double *X, int incX, int offset, int chunksize)
{

    // Get pointer to vector partition
    double* xVectorPart = &X[offset];

    // calculate dscal for determined part of vector X
    cblas_dscal(chunksize, Alpha, xVectorPart, 1);

    // We gather all parts of X and write it to X and get the result
    MPI_Allgather(xVectorPart, chunksize, MPI_DOUBLE, X, chunksize, MPI_DOUBLE, MPI_COMM_WORLD);

}

// Parallel cblas daxpy
// performs daxpy operation in parallel
inline void parallel_cblas_daxpy(int Size, double Alpha, double* X, int incX, double* Y, int incY, int offset, int chunksize)
{
    // Get pointer to vector partitions
    double* xVectorPart = &X[offset];
    double* yVectorPart = &Y[offset];

    // calculate daxpy for determined parts of X and Y vectors
    cblas_daxpy(chunksize, Alpha, xVectorPart, 1, yVectorPart, 1);

    // We gather all parts of X and write it to Y and get the result
    MPI_Allgather(yVectorPart, chunksize, MPI_DOUBLE, Y, chunksize, MPI_DOUBLE, MPI_COMM_WORLD);

}

inline void forward_sweep_parallel(double *A, int mSize, double *b, double *x, int offset, int size, int chunksize)
{
    // set x values to 0
    std::memset(x, 0, sizeof(double)*mSize);
    /*
    for(int i=0; i<mSize; i++)
    {
        for(int j=i-1; j>=0; j--)
            x[i] += A[i * mSize + j]*x[j];

        x[i] = (b[i] - x[i])/A[i* mSize + i];
    } */
}


/**
 * Parallel version of gauss_seidel_seq function
 */
inline int gauss_seidel_parallel(double *A, double *b, double *x0, size_t size, size_t maxIterations, double tolerance, int offset, int chunksize)
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

        // This cblas call calculates (b - U*x) and writes the result to b_new in parallel
        parallel_cblas_dgemv(CBLAS_ORDER::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans, size, size, -1.0, U, size, x, 1, 1.0, b_new, 1, offset, chunksize, size);

        // get new x vector
        forward_sweep_seq(L, size, b_new, x_new);

        // increase iteration count
        iterations++;
        
        // set elements of x_new to diff vector
        std::memcpy(diff, x_new, sizeof(double)*size);

        // get difference vector between x and x_new
        // and write it to diff vector
        parallel_cblas_daxpy(size, -1, x, 1, diff, 1, offset, chunksize);

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


inline int chebyshev_parallel(double *A, double *b, double *x0, size_t size, int maxIterations, double tolerance, double lMax, double lMin, int offset, int chunksize)
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
    parallel_cblas_dgemv(CBLAS_ORDER::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans, size, size, -1.0, A, size, x, 1, 1.0, r, 1, offset, chunksize, size);    

    // this keeps track of how many iterations have passed
    int iterations = 0;

    for(int i=0; i<maxIterations; i++)
    {
        iterations++;

        // fill x_pre with zeros
        std::memcpy(x_preNew,x_pre,sizeof(double)*size);
        gauss_seidel_parallel(preCond, r, x_preNew, size, 100, 1e-15, offset, chunksize);

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
            parallel_cblas_daxpy(size, beta, p, 1, p, 1, offset, chunksize);            
        }
        else
        {
            beta  = std::pow(c*alpha/2,2);
            alpha = 1/(d - beta/alpha);

            // These cblas calls calculate (p = z + beta*p) and writes the result to p
            parallel_cblas_dscal(size, beta, p , 1, offset, chunksize);
            parallel_cblas_daxpy(size, 1, x_preNew, 1, p, 1, offset, chunksize);
        }

        // x = x + alpha * p
        parallel_cblas_daxpy(size, alpha, p, 1, x, 1, offset, chunksize);

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
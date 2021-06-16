#include <iostream>
#include <cblas.h>
#include <utils.h>
#include <mpi.h>
#include <RootDir.h>
#include <string>


#include <matio.h>


#define tS(x) std::cout<<"\t"<<(#x)<<" == "<<(x)<<"\n"

#define COLS  12
#define ROWS  8

/*
int main(int argc, char** argv)
{

   enum CBLAS_ORDER order;
   enum CBLAS_TRANSPOSE transa;

   double *x, *y;

   double alpha, beta;
   int m, n, lda, incx, incy, i;

   order = CblasColMajor;
   transa = CblasNoTrans;

   m = 4; // Size of Column ( the number of rows ) 
   n = 4; // Size of Row ( the number of columns )
   lda = 4; // Leading dimension of 5 * 4 matrix is 5 
   incx = 1;
   incy = 1;
   alpha = 1;
   beta = 0;

   double a[m][n];

   //a = (double *)malloc(sizeof(double)*m*n);
   x = (double *)malloc(sizeof(double)*n);
   y = (double *)malloc(sizeof(double)*n);
   // The elements of the first column 
   a[0][0] = 1;
   a[0][1] = 2;
   a[0][2] = 3;
   a[0][3] = 4;
   // The elements of the second column 
   a[1][0] = 1;
   a[1][1] = 1;
   a[1][2] = 1;
   a[1][3] = 1;
   // The elements of the third column
   a[2][0] = 3;
   a[2][1] = 4;  
   a[2][2] = 5;
   a[2][3] = 6;
   // The elements of the fourth column
   a[3][0] = 5;
   a[3][1] = 6;
   a[3][2] = 7;
   a[3][3] = 8;
   // The elemetns of x and y
   x[0] = 1;
   x[1] = 2;
   x[2] = 1;
   x[3] = 1;
   y[0] = 0;
   y[1] = 0;
   y[2] = 0;
   y[3] = 0;
   
   cblas_dgemv( order, transa, m, n, alpha, *a, lda, x, incx, beta,
                y, incy );
   cblas_dscal(n, 2, y, 1);
   // Print y
   //for( i = 0; i < n; i++ ) 
   //   printf(" y%d = %f\n", i, y[i]);
   //free(a);

   //printMat(*a, m, n, MatOrder::colMajor);
   printVec(y, n);
   free(x);
   free(y);
   return 1;

}*/


int main(int argc, char **argv) {

    

    std::string filename = ROOT_DIR + std::string("matrices/3x3IdentityMatrix.mat");

    mat_t *mat = Mat_Open(filename.c_str(), MAT_ACC_RDONLY);

    matvar_t *matVar = Mat_VarRead(mat, (char*)"A");

            unsigned xSize = matVar->nbytes/matVar->data_size ;
            const double *xData = static_cast<const double*>(matVar->data) ;
            for(int i=0; i<xSize; ++i)
            {
                std::cout<<"\tx["<<i<<"] = "<<xData[i]<<"\n" ;
            }
            std::cout<<"\n" ;
            for(int i=0; i<matVar->rank; ++i)
            {
                std::cout<<"\tdim["<<i<<"] == "<<matVar->dims[i]<<"\n" ;
            }
    

    // MPI_Init(&argc, &argv);
    // int p, rank;
    // MPI_Comm_size(MPI_COMM_WORLD, &p);
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // char i;

    // if(rank == 0)
    // {
    //     std::cout << "Hello Beginning" << std::endl;
    // }

    // char a[ROWS*COLS];
    // const int NPROWS=2;  /* number of rows in _decomposition_ */
    // const int NPCOLS=3;  /* number of cols in _decomposition_ */
    // const int BLOCKROWS = ROWS/NPROWS;  /* number of rows in _block_ */
    // const int BLOCKCOLS = COLS/NPCOLS; /* number of cols in _block_ */

    // if (rank == 0) {
    //     for (int ii=0; ii<ROWS*COLS; ii++) {
    //         a[ii] = (char)ii;
    //     }
    // }

    // if (p != NPROWS*NPCOLS) {
    //     fprintf(stderr,"Error: number of PEs %d != %d x %d\n", p, NPROWS, NPCOLS);
    //     MPI_Finalize();
    //     exit(-1);
    // }
    // char b[BLOCKROWS*BLOCKCOLS];
    // for (int ii=0; ii<BLOCKROWS*BLOCKCOLS; ii++) b[ii] = 0;

    // MPI_Datatype blocktype;
    // MPI_Datatype blocktype2;

    // MPI_Type_vector(BLOCKROWS, BLOCKCOLS, COLS, MPI_CHAR, &blocktype2);
    // MPI_Type_create_resized( blocktype2, 0, sizeof(char), &blocktype);
    // MPI_Type_commit(&blocktype);

    // int disps[NPROWS*NPCOLS];
    // int counts[NPROWS*NPCOLS];
    // for (int ii=0; ii<NPROWS; ii++) {
    //     for (int jj=0; jj<NPCOLS; jj++) {
    //         disps[ii*NPCOLS+jj] = ii*COLS*BLOCKROWS+jj*BLOCKCOLS;
    //         counts [ii*NPCOLS+jj] = 1;
    //     }
    // }

    // MPI_Scatterv(a, counts, disps, blocktype, b, BLOCKROWS*BLOCKCOLS, MPI_CHAR, 0, MPI_COMM_WORLD);
    // /* each proc prints it's "b" out, in order */
    // for (int proc=0; proc<p; proc++) {
    //     if (proc == rank) {
    //         printf("Rank = %d\n", rank);
    //         if (rank == 0) {
    //             printf("Global matrix: \n");
    //             for (int ii=0; ii<ROWS; ii++) {
    //                 for (int jj=0; jj<COLS; jj++) {
    //                     printf("%3d ",(int)a[ii*COLS+jj]);
    //                 }
    //                 printf("\n");
    //             }
    //         }
    //         printf("Local Matrix:\n");
    //         for (int ii=0; ii<BLOCKROWS; ii++) {
    //             for (int jj=0; jj<BLOCKCOLS; jj++) {
    //                 printf("%3d ",(int)b[ii*BLOCKCOLS+jj]);
    //             }
    //             printf("\n");
    //         }
    //         printf("\n");
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    // MPI_Finalize();

    return 0;
}
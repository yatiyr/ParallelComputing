#include <iostream>
#include <cblas.h>
#include <Utils.h>
#include <RootDir.h>
#include <string>
#include <ParallelFunctions.h>
#include <SequentialFunctions.h>
#include <fstream>
#include <sstream>
#include <MatfileReader.h>
#include <Timer.h>
#include <matio.h>


#include <mpi.h>

#define NOT_ENOUGH_PROCESSES_NUM_ERROR 1
#define ROOT_PROCESS 0

int main(int argc, char **argv) {

    // each process will know its
    // number and process size from
    // these variables
    int worldSize;
    int worldRank;

    // we also need an offset value
    // for partitioning matrices and
    // vectors
    int offset;

    // we need to determine size of
    // data chunks as well
    int chunkSize;

    // Root process will read matrix and
    // fill in this pointer
    // and broadcast it to other processes
    // with its size of course
    // we also have method name, method namesize
    // and max,min eigenvalues of our matrix
    double* Matrix = nullptr;
    int MatrixDim;
    char* MethodName = nullptr;
    int methodNameSize;
    double maxEig, minEig;


    MPI_Init(&argc, &argv);

    // Each process gets how many processors are working
    // and their number
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    // we need more than 2 and even number of processes
    if(worldSize < 2 && worldSize % 2 == 0)
    {
        MPI_Abort(MPI_COMM_WORLD, NOT_ENOUGH_PROCESSES_NUM_ERROR);        
    }



    // Root process reads config file and forms the
    // matrix we are going to work with. 
    // 
    // Then, it broadcasts size of the matrix and
    // the data
    if(worldRank == ROOT_PROCESS)
    {
        // Read config file
        // TODO: CAN MOVE THIS CODE INTO A CLASS OR MAKE IT A FUNCTION I DON'T KNOW :D
        std::string configPath = std::string(ROOT_DIR) + "config/config.conf";
        std::ifstream file(configPath);
        std::string line = "";


        std::string matrixName;
        std::string rightHandName;
        std::string methodName;
        double lMax, lMin;

        while(std::getline(file,line))
        {
            std::stringstream ss(line);
            std::string word = "";
            ss >> word;
            
            if(word == "matrix")
            {
                ss >> word;
                matrixName = word;
                ss >> lMax >> lMin;
            }
            else if(word == "rightHand")
            {
                ss >> word;
                rightHandName = word;
            }
            else if(word == "it_method")
            {
                ss >> word;
                methodName = word;
            }
        }

        std::string matrixPath = std::string(ROOT_DIR) + "matfiles/" + matrixName;
        std::string rightHandPath = std::string(ROOT_DIR) + "matfiles/" + rightHandName;

        // Fill in work matrix and its dimension
        Mat matrix = MatfileReader::ReadMat(matrixName);
        Matrix      = matrix.data;
        MatrixDim   = matrix.columnSize;

        MethodName = (char*)malloc(sizeof(char)*(methodName.size() + 1));
        std::copy(methodName.begin(), methodName.end(), MethodName);
        MethodName[methodName.size()] = '\0';

        methodNameSize = methodName.size();

        maxEig = lMax;
        minEig = lMin;

        // Broadcast this data to other processes so they can read and
        // use it too
        MPI_Bcast(&MatrixDim, 1, MPI_INT, ROOT_PROCESS, MPI_COMM_WORLD);
        MPI_Bcast(Matrix, MatrixDim*MatrixDim, MPI_DOUBLE, ROOT_PROCESS, MPI_COMM_WORLD);
        MPI_Bcast(&methodNameSize, 1, MPI_INT, ROOT_PROCESS, MPI_COMM_WORLD);        
        MPI_Bcast(MethodName, methodName.size() + 1, MPI_CHAR, ROOT_PROCESS, MPI_COMM_WORLD);
        MPI_Bcast(&maxEig, 1, MPI_DOUBLE, ROOT_PROCESS, MPI_COMM_WORLD);
        MPI_Bcast(&minEig, 1, MPI_DOUBLE, ROOT_PROCESS, MPI_COMM_WORLD);

    }
    else
    {
        // We use bcast in other processes as well in order to get the data
        // we first get the dimensions
        MPI_Bcast(&MatrixDim, 1, MPI_INT, ROOT_PROCESS, MPI_COMM_WORLD);

        // after we get size of matrix, we allocate space for our Matrix pointer
        Matrix = (double*)malloc(sizeof(double)*MatrixDim*MatrixDim);
        MPI_Bcast(Matrix, MatrixDim*MatrixDim, MPI_DOUBLE, ROOT_PROCESS, MPI_COMM_WORLD);

        // Get method name but first allocate
        MPI_Bcast(&methodNameSize, 1, MPI_INT, ROOT_PROCESS, MPI_COMM_WORLD);        
        MethodName = (char*)malloc(sizeof(char)*(methodNameSize + 1));
        MPI_Bcast(MethodName, methodNameSize + 1, MPI_CHAR, ROOT_PROCESS, MPI_COMM_WORLD);

        // Get max min eigenvalues
        MPI_Bcast(&maxEig, 1, MPI_DOUBLE, ROOT_PROCESS, MPI_COMM_WORLD);
        MPI_Bcast(&minEig, 1, MPI_DOUBLE, ROOT_PROCESS, MPI_COMM_WORLD);        

    }

    
    // Allocate initial guess and right hand side vectors
    double* x0 = (double*) calloc(MatrixDim, sizeof(double));
    double* b  = ones(MatrixDim);

    // calculate chunk size and offset
    chunkSize = MatrixDim / worldSize;
    offset = worldRank * chunkSize;

    int iterations = 0;

    {
        if(std::strcmp(MethodName, "gauss-seidel") == 0)
        {
            iterations = gauss_seidel_parallel(Matrix, b, x0, MatrixDim, 100, 1.5e-10, offset, chunkSize);
        }
        else if(std::strcmp(MethodName, "jacobi") == 0)
        {
            iterations = jacobi_seq(Matrix, b, x0, MatrixDim, 100, 1.5e-10);
        }
        else if(std::strcmp(MethodName, "chebyshev") == 0)
        {
            iterations = chebyshev_parallel(Matrix, b, x0, MatrixDim, 100, 1.5e-10, maxEig, minEig, offset, chunkSize);
        }
        else
        {
            printf("%s\n",MethodName);
        }

    }

    if(worldRank == ROOT_PROCESS)
        std::cout << "Root Process reports that " << iterations << " iterations have passed." << '\n';

    MPI_Finalize();


//    double* x0 = (double*) calloc(matrix.columnSize, sizeof(double));
//    double* b  = ones(matrix.columnSize);




    /*

    int iterations = 0;

    {
        Timer t;
        if(methodName == "gauss-seidel")
        {
            iterations = gauss_seidel_seq(matrix.data, b, x0, matrix.rowSize, 100, 1.5e-10);
        }
        else if(methodName == "jacobi")
        {
            iterations = jacobi_seq(matrix.data, b, x0, matrix.rowSize, 100, 1.5e-10);
        }
        else if(methodName == "chebyshev")
        {
            iterations = chebyshev_seq(matrix.data, b, x0, matrix.rowSize, 100, 1.5e-10, lMax, lMin);
        }

    }

    std::cout << iterations << " iterations have passed." << '\n';


    free(x0);
    free(b);

    */








    /*

    MPI_Init(nullptr, nullptr);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    printf("Hello world from processor %s, rank %d out of %d processors \n", processor_name, world_rank, world_size);

    MPI_Finalize(); */


    // std::string filename = ROOT_DIR + std::string("matrices/1138_bus.mat");

    // mat_t *mat = Mat_Open(filename.c_str(), MAT_ACC_RDWR);

    // matvar_t *matVar = Mat_VarRead(mat, (char*)"Problem");

    //         unsigned xSize = matVar->nbytes/matVar->data_size ;
    //         matvar_t *xData = matVar->data;

            /*
            for(int i=0; i<xSize; ++i)
            {
                std::cout<<"\tx["<<i<<"] = "<<xData[i]<<"\n" ;
            }
            std::cout<<"\n" ;
            for(int i=0; i<matVar->rank; ++i)
            {
                std::cout<<"\tdim["<<i<<"] == "<<matVar->dims[i]<<"\n" ;
            } */
    

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
#ifndef __MATFILE_READER_H__
#define __MATFILE_READER_H__

#include <matio.h>
#include <RootDir.h>
#include <string>
#include <Utils.h>


struct Mat
{
    int rowSize;
    int columnSize;
    double* data;
};

class MatfileReader
{
public:

    static Mat ReadMat(std::string matName)
    {
        Mat result;

        std::string path = std::string(ROOT_DIR) + "matfiles/" + matName;
        std::string varName = matName.substr(0, matName.find("."));

        mat_t *mat = Mat_Open(path.c_str(), MAT_ACC_RDWR);
        matvar_t *matVar = Mat_VarRead(mat, varName.c_str());

        result.rowSize    = (int) matVar->dims[0];
        result.columnSize = (int) matVar->dims[1];
        result.data = (double*) matVar->data;

        // Matio reads matrices column major order, I wanted to change it
        // to row major order
        changeCol2Row(result.data, result.rowSize, result.columnSize);
        
        return result;
    }
};

#endif
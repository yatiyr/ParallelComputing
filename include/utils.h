#ifndef __UTILS_H__
#define __UTILS_H__

#include <iostream>

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


#endif
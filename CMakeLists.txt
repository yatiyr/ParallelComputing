cmake_minimum_required(VERSION 3.10)
project(ParallelComputing VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")



file(GLOB_RECURSE SOURCE_FILES 
     ${CMAKE_SOURCE_DIR}/src/*.c
     ${CMAKE_SOURCE_DIR}/src/*.cpp)

file(GLOB_RECURSE HEADER_FILES
     ${CMAKE_SOURCE_DIR}/src/*.h
     ${CMAKE_SOURCE_DIR}/src/*.hpp)

configure_file(src/helpers/RootDir.h.in src/RootDir.h)
include_directories(${CMAKE_BINARY_DIR}/src)

add_executable(${PROJECT_NAME} ${HEADER_FILES} ${SOURCE_FILES})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(CBLAS REQUIRED)
message(STATUS "CBLAS included at ${CBLAS_INCLUDE_DIR}")

find_package(BLAS REQUIRED)

find_package(MPI REQUIRED)

find_package(MATIO REQUIRED)

include_directories(
    "${CMAKE_SOURCE_DIR}/src"
    "${CMAKE_SOURCE_DIR}/include"
    "{CMAKE_SOURCE_DIR}/thirdparty/include"
    ${CBLAS_INCLUDE_DIR}
    ${MATIO_INCLUDE_DIR}
    MPI::MPI_CXX_INCLUDE_PATH
)

message(STATUS "BLAS Library at ${BLAS_LIBRARY}")

set(LIBS pthread ${CBLAS_LIBRARY} ${BLAS_LIBRARY} ${MATIO_LIBRARY} gfortran MPI::MPI_CXX)

target_link_libraries(${PROJECT_NAME} ${LIBS})
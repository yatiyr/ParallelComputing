# - Try to find CBLAS
# Once done, this will define
#
# CBLAS_FOUND - system has CBLAS
# CBLAS_INCLUDE_DIR - the CBLAS include directories
# CBLAS_LIBRARIES - link these to use CBLAS
FIND_PATH(CBLAS_INCLUDE_DIR cblas.h
	/usr/include
	/usr/local/include
	/opt/local/include
    "${CMAKE_SOURCE_DIR}/include"
    "${CMAKE_SOURCE_DIR}/thirdparty/include"
)
FIND_LIBRARY( CBLAS_LIBRARY cblas
	/usr/lib64
	/usr/lib
	/usr/local/lib
	/opt/local/lib
    "${CMAKE_SOURCE_DIR}/lib"
    "${CMAKE_SOURCE_DIR}/thirdparty/lib"
)
IF(CBLAS_INCLUDE_DIR AND CBLAS_LIBRARY)
	SET( CBLAS_FOUND TRUE )
	SET( CBLAS_LIBRARIES ${CBLAS_LIBRARY} )
ENDIF(CBLAS_INCLUDE_DIR AND CBLAS_LIBRARY)
IF(CBLAS_FOUND)
	IF(NOT CBLAS_FIND_QUIETLY)
	MESSAGE(STATUS "Found CBLAS: ${CBLAS_LIBRARY}")
	ENDIF(NOT CBLAS_FIND_QUIETLY)
ELSE(CBLAS_FOUND)
	IF(CBLAS_FIND_REQUIRED)
	MESSAGE(FATAL_ERROR "Could not find libcblas")
	ENDIF(CBLAS_FIND_REQUIRED)
ENDIF(CBLAS_FOUND)
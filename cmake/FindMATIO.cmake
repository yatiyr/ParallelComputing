# - Try to find MATIO
# Once done, this will define
#
# MATIO_FOUND - system has MATIO
# MATIO_INCLUDE_DIR - the MATIO include directories
# MATIO_LIBRARIES - link these to use MATIP
FIND_PATH(MATIO_INCLUDE_DIR matio.h
	/usr/include
	/usr/local/include
	/opt/local/include
    "${CMAKE_SOURCE_DIR}/include"
    "${CMAKE_SOURCE_DIR}/thirdparty/include"
)
FIND_LIBRARY(MATIO_LIBRARY matio
	/usr/lib64
	/usr/lib
	/usr/local/lib
	/opt/local/lib
    "${CMAKE_SOURCE_DIR}/lib"
    "${CMAKE_SOURCE_DIR}/thirdparty/lib"
)
IF(MATIO_INCLUDE_DIR AND MATIO_LIBRARY)
	SET( MATIO_FOUND TRUE )
	SET( MATIO_LIBRARIES ${MATIO_LIBRARY} )
ENDIF(MATIO_INCLUDE_DIR AND MATIO_LIBRARY)
IF(MATIO_FOUND)
	IF(NOT MATIO_FIND_QUIETLY)
	MESSAGE(STATUS "Found MATIO: ${MATIO_LIBRARY}")
	ENDIF(NOT MATIO_FIND_QUIETLY)
ELSE(MATIO_FOUND)
	IF(MATIO_FIND_REQUIRED)
	MESSAGE(FATAL_ERROR "Could not find libmatio")
	ENDIF(MATIO_FIND_REQUIRED)
ENDIF(MATIO_FOUND)
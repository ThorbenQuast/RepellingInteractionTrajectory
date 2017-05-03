#nice tutorial: http://www.cmake.org/cmake-tutorial/

cmake_minimum_required (VERSION 2.6)
project (StrongInteractionTrajectory)
# The version number.
set (StrongInteractionTrajectory_VERSION_MAJOR 0)
set (StrongInteractionTrajectory_VERSION_MINOR 1)
 
# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/config/StrongInteractionTrajectoryConfig.in"
  "${PROJECT_BINARY_DIR}/config/StrongInteractionTrajectoryConfig.h"
  )
 

### LIBRARY AND TOOLS #######################################################
# Add PXL and link it
include_directories (
	${PXL_SOURCE_DIR}/modules/include
	${PXL_SOURCE_DIR}/scripting/include
	${PXL_SOURCE_DIR}/algorithms/include
	${PXL_SOURCE_DIR}/core/include
	${PXL_SOURCE_DIR}/hep/include
	${PXL_SOURCE_DIR}/astro/include
	${PXL_SOURCE_DIR}/xml/include
	${PXL_SOURCE_DIR}/lib
	${CMAKE_BINARY_DIR}/include
)

# library name
SET( libname ${PROJECT_NAME} )

# include directories
INCLUDE_DIRECTORIES( ./include )
INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}/config")




# add library the other includes as directory
AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_LIBRARY( ${libname} ${library_sources} )
target_link_libraries (${libname} pxl-core pxl-hep)

LINK_LIBRARIES( ${libname} )

# make the executable main
ADD_EXECUTABLE( trajectory trajectory.cc )
ADD_EXECUTABLE( integral integral.cc )







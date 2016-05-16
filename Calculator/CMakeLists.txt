cmake_minimum_required(VERSION 2.8)
 
PROJECT(FinalProject)

set(CMAKE_BUILD_TYPE Release)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_library(RegionGrowingNoThreshold lib/RegionGrowingNoThreshold.cxx lib/RegionGrowingNoThreshold.h)

add_executable(FinalProject MACOSX_BUNDLE FinalProject)

add_definitions(-std=c++11)

if(VTK_LIBRARIES)
  target_link_libraries(FinalProject ${VTK_LIBRARIES})
else()
  target_link_libraries(FinalProject vtkHybrid vtkWidgets)
endif()

#Root from CERN
#Needed for advance curvature fittin, not used
#list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#find_package(ROOT REQUIRED COMPONENTS RIO Net)
#include(${ROOT_USE_FILE})
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
#include_directories(${ROOT_INCLUDE_DIRS})

file( GLOB LIB_SOURCES lib/*.cpp )
file( GLOB LIB_HEADERS lib/*.h )
add_library( YourLib ${LIB_SOURCES} ${LIB_HEADERS} )

find_package(ITK REQUIRED)
FIND_PACKAGE(GSL REQUIRED)

include(${ITK_USE_FILE})

#if (ITKVtkGlue_LOADED)
#  find_package(VTK REQUIRED)
#  include(${VTK_USE_FILE})
#else()
#  find_package(ItkVtkGlue REQUIRED)
#  include(${ItkVtkGlue_USE_FILE})
#  set(Glue ItkVtkGlue)
#endif()
 
target_link_libraries( FinalProject YourLib RegionGrowingNoThreshold ${VTK_LIBRARIES} ${ITK_LIBRARIES} ${ROOT_LIBRARIES} ${GSL_LIBRARIES})

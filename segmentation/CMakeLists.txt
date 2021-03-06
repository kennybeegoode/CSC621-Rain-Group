set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
cmake_minimum_required(VERSION 2.8)

project(ConnectedThresholdImageFilter)
project(ConnectedThresholdImageFilterNoThreshold)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
if (ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
else()
  find_package(ItkVtkGlue REQUIRED)
  include(${ItkVtkGlue_USE_FILE})
  set(Glue ItkVtkGlue)
endif()

add_executable(ConnectedThresholdImageFilter MACOSX_BUNDLE ConnectedThresholdImageFilter.cxx)
target_link_libraries(ConnectedThresholdImageFilter
  ${Glue}  ${VTK_LIBRARIES} ${ITK_LIBRARIES})

add_executable(ConnectedThresholdImageFilterNoThreshold MACOSX_BUNDLE ConnectedThresholdImageFilterNoThreshold.cxx)
  target_link_libraries(ConnectedThresholdImageFilterNoThreshold
    ${Glue}  ${VTK_LIBRARIES} ${ITK_LIBRARIES})

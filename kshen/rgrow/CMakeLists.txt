cmake_minimum_required(VERSION 2.8)

project(rgrow)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(rgrow rgrow.cxx )

target_link_libraries(rgrow ${ITK_LIBRARIES})
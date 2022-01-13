cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(brufit)

	
set(CMAKE_INSTALL_BINDIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(CMAKE_CXX_FLAGS "-fPIC -O3")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}) ############USEROOTSYS

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS Proof RooStats MathMore)  ###########USEROOTSYS

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE}) ##########USEROOTSYS
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/core)
add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/core)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/coreclang)
#add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/coreclang)


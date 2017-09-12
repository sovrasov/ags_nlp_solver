cmake_minimum_required(VERSION 3.1)

set(PROJECT_NAME_STR ags)
project(${PROJECT_NAME_STR})

#initial setup
enable_language(C)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)

#build optimizer
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/core/include)

file(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/core/src/*.cpp
${CMAKE_CURRENT_SOURCE_DIR}/core/include/*.hpp)

add_library(${PROJECT_NAME_STR} ${SRC_FILES})
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  target_compile_options(${PROJECT_NAME_STR} INTERFACE -Wall -fopenmp -march=native)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  add_definitions(/EHsc /openmp)
endif()

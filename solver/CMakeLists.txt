cmake_minimum_required(VERSION 3.5)

set(PROJECT_NAME_STR nlp_solver)
project(${PROJECT_NAME_STR})

enable_language(C)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)

file(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
${CMAKE_CURRENT_SOURCE_DIR}/core/include/*.hpp)

add_library(${PROJECT_NAME_STR} ${SRC_FILES})

set_property(TARGET ${PROJECT_NAME_STR} PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(${PROJECT_NAME_STR} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${GO_INTERFACE_INCLUDE_DIR}/)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  target_compile_options(${PROJECT_NAME_STR} INTERFACE -Wall -fopenmp -march=native)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  add_definitions(/EHsc /openmp)
endif()

if (${BUILD_BINDINGS})
  pybind11_add_module(ags_solver ${CMAKE_CURRENT_SOURCE_DIR}/bindings_src/py_bindings.cpp)
  target_include_directories(ags_solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${GO_INTERFACE_INCLUDE_DIR}/)
  target_link_libraries(ags_solver PUBLIC ${PROJECT_NAME_STR})
endif()

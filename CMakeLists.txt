cmake_minimum_required(VERSION 3.5)

add_subdirectory(3rd-party/global-optimization-test-problems)
set(CMDLINE_INCLUDE ${CMAKE_CURRENT_LIST_DIR}/3rd-party/cmdline)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

option(BUILD_BINDINGS "Build python bindings" OFF)
if (${BUILD_BINDINGS})
  add_subdirectory(3rd-party/pybind11)
endif()

add_subdirectory(solver)
add_subdirectory(samples)

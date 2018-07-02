cmake_minimum_required(VERSION 3.5)

add_subdirectory(solver)

set(CMDLINE_INCLUDE "3rd-party/cmdline")
add_subdirectory(3rd-party/global-optimization-test-problems)
add_subdirectory(sample)

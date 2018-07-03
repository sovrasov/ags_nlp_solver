cmake_minimum_required(VERSION 3.5)

add_subdirectory(3rd-party/global-optimization-test-problems)
set(CMDLINE_INCLUDE ${CMAKE_CURRENT_LIST_DIR}/3rd-party/cmdline)

add_subdirectory(solver)
add_subdirectory(sample)

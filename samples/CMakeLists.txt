cmake_minimum_required(VERSION 3.5)

set(PROJECT_NAME_STR sample)
project(${PROJECT_NAME_STR})

enable_language(C)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)

add_executable(solve_set ${CMAKE_CURRENT_SOURCE_DIR}/solve_set.cpp)
target_include_directories(solve_set PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMDLINE_INCLUDE})
target_link_libraries(solve_set nlp_solver gkls_functions grishagin_functions)

add_executable(solve_constrained ${CMAKE_CURRENT_SOURCE_DIR}/constrained.cpp)
target_link_libraries(solve_constrained nlp_solver gkls_functions grishagin_functions)

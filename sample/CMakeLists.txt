cmake_minimum_required(VERSION 3.5)

set(PROJECT_NAME_STR sample)
project(${PROJECT_NAME_STR})

enable_language(C)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)

file(GLOB SAMPLE_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp ${CMAKE_CURRENT_SOURCE_DIR}/*.h)

add_executable(${PROJECT_NAME_STR} ${SAMPLE_SRC_FILES})
target_include_directories(${PROJECT_NAME_STR} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(${PROJECT_NAME_STR} nlp_solver)
target_link_libraries(${PROJECT_NAME_STR} gkls_functions)

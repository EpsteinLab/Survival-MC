cmake_minimum_required(VERSION 3.8)
project(MC_simulation)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

include_directories("${PROJECT_SOURCE_DIR}")

add_subdirectory(src "${PROJECT_SOURCE_DIR}/src")

add_executable(MC_simulation "${PROJECT_SOURCE_DIR}/main.cpp")

target_link_libraries(MC_simulation PUBLIC src)

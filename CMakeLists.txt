cmake_minimum_required(VERSION 3.15)

project(mse)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3 -g -Wall")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
aux_source_directory(. DIR_SRCS)
add_executable(mse ${DIR_SRCS})

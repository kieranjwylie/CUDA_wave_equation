cmake_minimum_required(VERSION 3.18)
project(wave LANGUAGES CXX CUDA)

set(CMAKE_CUDA_STANDARD 14)

file(GLOB SOURCES
     src/*.cpp
     src/*.cu
     src/*.hpp
     src/*.cuh)

add_executable(wave ${SOURCES})

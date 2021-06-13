cmake_minimum_required(VERSION 3.10)
project(inv)

set(CMAKE_CXX_STANDARD 14)

find_package(MPI REQUIRED)
include_directories("${MPI_INCLUDE_PATH}")

add_executable(inv main.cpp
               adjugate.cpp
               adjugate.h
               elimination.cpp
               elimination.h)

target_link_libraries(inv ${MPI_CXX_LIBRARIES})
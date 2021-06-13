cmake_minimum_required(VERSION 3.10)
project(inv)

set(CMAKE_CXX_STANDARD 14)

add_executable(inv main.cpp
               adjugate.cpp
               adjugate.h
               elimination.cpp
               elimination.h)
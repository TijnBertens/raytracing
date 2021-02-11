cmake_minimum_required(VERSION 3.12)
project(raytracing)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)

add_executable(raytracing src/main.cpp)
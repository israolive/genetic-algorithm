cmake_minimum_required(VERSION 3.16)
project(genetic_algorithm)

set(CMAKE_CXX_STANDARD 17)

include_directories(.)

add_executable(genetic_algorithm
        model/genetic_algorithm.hpp
        model/stochastic.hpp
        main.cpp)

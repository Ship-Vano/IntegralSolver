cmake_minimum_required(VERSION 3.27)
project(Code)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(InputData)
include_directories(Libs)
include_directories(OutputData)

add_executable(Code
        Libs/algebra.h
        Libs/FileIO.h
        main.cpp
        Libs/IntegralProblem.h
        Libs/IntegralProblemSolver.h
        Libs/IntegralProblem.cpp
        Libs/IntegralProblemSolver.cpp
        Libs/algebra.cpp
        Libs/SLAEsolver.h)

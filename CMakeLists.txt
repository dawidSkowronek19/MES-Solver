cmake_minimum_required(VERSION 3.15)
project(FemSolver1D VERSION 1.0 LANGUAGES CXX)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


find_package(GSL REQUIRED)

find_package(Eigen3 REQUIRED) 


add_executable(solver_exe 
    linear.cpp 
    Solver.cpp 
    SolverPhysics.cpp 
    MathSolver.cpp
    Grid.cpp
    ShapeFunction.cpp
)


target_link_libraries(solver_exe PRIVATE GSL::gsl GSL::gslcblas Eigen3::Eigen)


target_precompile_headers(solver_exe PRIVATE <Eigen/Dense> <vector> <iostream> <cmath>)

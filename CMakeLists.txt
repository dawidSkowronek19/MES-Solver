cmake_minimum_required(VERSION 3.15)
project(FemSolver1D VERSION 1.0 LANGUAGES CXX)

# Ustawienie standardu C++ (np. C++17 dla strukturyzowanych powiązań auto [A, B])
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 1. Znajdź wymagane biblioteki zewnętrzne
find_package(GSL REQUIRED)
# Eigen często dostarcza konfigurację dla CMake pod nazwą Eigen3
find_package(Eigen3 REQUIRED) 

# 2. Zdefiniuj plik wykonywalny i wskaż pliki źródłowe
add_executable(solver_exe 
    linear.cpp 
    Solver.cpp 
    SolverPhysics.cpp 
    MathSolver.cpp
    Grid.cpp
    ShapeFunction.cpp
)

# 3. Podepnij biblioteki i ścieżki do plików nagłówkowych pod nasz program
target_link_libraries(solver_exe PRIVATE GSL::gsl GSL::gslcblas Eigen3::Eigen)

# 4. Szybka kompilacja: Dodaj przedkompilowane nagłówki dla Eigena i STL
target_precompile_headers(solver_exe PRIVATE <Eigen/Dense> <vector> <iostream> <cmath>)

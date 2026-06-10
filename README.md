# C++ 1D Finite Element Method (FEM) Solver

This repository hosts a highly modular and robust 1D Finite Element Method (FEM) solver written in C++. The ultimate goal of this project is to build a highly optimized, user-friendly tool for solving stationary, non-linear, time-dependent partial differential equations (PDEs), and eigenvalue problems.

## Current Status
The current stable baseline version is in branch `main`.

## Compile Instructions
To compile the project using CMake, run the following command in the root directory:
```bash
cmake --build ./build/ -j"X" 
```
*(Where `X` = the number of CPU cores used for compiling)*

> **Note:** This project utilizes the **Eigen3** library for algebraic operations and the **GNU Scientific Library (GSL)** for numerical integration. Ensure both are installed on your system before building.

---

## What's New in Version 0.2

**1. Full Solver Capabilities Expanded**
The solver engine has been vastly expanded to handle a wide variety of physical equations:
* **Time-Dependent Solver:** A new solver utilizing the Newmark-beta method for integrating 1D time-dependent differential equations.
* **Eigenvalue Solver:** Implemented capabilities to solve both General and Symmetric Eigenvalue problems using Eigen's `GeneralizedSelfAdjointEigenSolver` and `GeneralizedEigenSolver`. (eigen.py requires to set which mode you want to plot)
* **Nonlinear Stationary Solver:** The Newton-Raphson iteration method is now fully implemented and operational, successfully resolving non-linearities with automated Tangent Matrix assembly.

**2. Higher-Order Shape Functions**
The implementation now supports arbitrary polynomial degrees for shape functions. The generic `ShapeFunction` class allows users to easily define the degree of the polynomials, enabling advanced p-FEM capabilities.

**3. Advanced Grid Generation**
The `Grid_1D` class has been heavily upgraded. It now supports non-uniform grid generation with specific attraction constants (`growFactor`, `A`) and utilizes Gaussian quadrature (via GSL) to compute precise node distributions based on boundary condition locations.

**4. Refined Architecture: MathSolver & Physics**
The strict Object-Oriented separation continues to be improved:
* **`MathSolver`:** Extended to handle LU decompositions (`PartialPivLU`) for general systems and Bunch-Kaufman (`LDLT`) for symmetric systems.
* **`Physics`:** An abstract class with specific implementations (`Poisson`, `GeneralPDE`, `GeneralSymetricPDE`) that fully encapsulate the mathematical models and their derivatives (Jacobians) for non-linear solving.

---

## Current Development Goals & Future Roadmap

Active development is focused on making the solver significantly faster and more accessible. The immediate roadmap includes:

* **User-Friendly API & Architecture:** Continuously refining the public API of the solver classes so that setting up physics, meshes, simulation parameters, and boundary conditions is as intuitive and clean as possible.
* **Input File Parser:** Developing a robust parser to read simulation parameters, boundary conditions, and physical constants directly from input files. This will eliminate the need to recompile the C++ source code for every new simulation setup.
* **OpenMP Parallelization:** Implementing multi-threading using OpenMP (`omp`). The focus will be on parallelizing the computationally heavy matrix and load vector assembly processes, as well as local integration loops, to drastically reduce computation times for large node counts.
* **General Code Optimization:** Profiling the codebase to reduce memory overhead and optimize memory access patterns during matrix population. Implementing sparse matrixes

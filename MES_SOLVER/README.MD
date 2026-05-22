# 1D FEM Solver - Baseline

This repository contains the **stable baseline version** of a 1D Finite Element Method solver. It includes fully working implementations for:
* **1D Stationary Analysis** (solving standard boundary value problems).
* **Modal Analysis** (solving generalized eigenvalue problems to find eigenvalues and eigenvectors).

## Prerequisites
To compile and run this code, you need a C++17 compatible compiler and the following libraries installed on your system:
* **Eigen3** (for matrix operations and algebra solvers)
* **GSL (GNU Scientific Library)** (for numerical integration/Gauss-Legendre quadrature)

**Ubuntu/Debian installation:**
```bash
sudo apt update
sudo apt install build-essential libgsl-dev libeigen3-dev


***********************************

**IMPORTANT**
.
├── main.cpp
├── Solver.cpp
├── Solver.h
├── Grid.cpp
├── Grid.h
└── output/       <-- YOU MUST CREATE THIS FOLDER

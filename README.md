# 1D FEM Solver

A C++ implementation of the 1D Finite Element Method for solving boundary value problems, transient dynamics, and eigenvalue problems. The solver reads configurations and mathematical expressions at runtime, eliminating the need to recompile when changing equation parameters.

## Capabilities

* **Stationary Linear:** Solves standard 1D linear PDEs.
* **Stationary Nonlinear:** Implements the Newton-Raphson method with exact analytical tangent matrix (Jacobian) assembly.
* **Time-Dependent:** Uses the Newmark-beta integration scheme for transient problems.
* **Eigenvalue Problems:** Computes natural frequencies and mode shapes via `Eigen::GeneralizedSelfAdjointEigenSolver`.
* **Runtime Parsing:** Uses `ExprTk` to parse continuous mathematical functions (e.g., source terms, initial conditions) directly from configuration files.
* **Adaptive Grid:** Supports uniform meshing and local mesh refinement using Gaussian distribution functions.

## Dependencies

* **C++17** compatible compiler (GCC, Clang, MSVC)
* **CMake** (>= 3.10)
* **Eigen3** (Linear algebra operations)
* **GSL - GNU Scientific Library** (Gauss-Legendre numerical integration)
* **Python 3** with `matplotlib` and `numpy` (For post-processing scripts)
* *Note: `exprtk.hpp` is included in the source tree.*

## Build Instructions

The project uses CMake. To build the executable from the source:

```bash
mkdir build && cd build
cmake ..
cmake --build . -j8

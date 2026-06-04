# C++ Finite Element Method (FEM) Solver

This repository hosts a Finite Element Method (FEM) solver written in C++. The ultimate goal of this project is to build a robust, highly optimized, and user-friendly tool for solving various differential equations.

## Current Status: Stable Baseline
The current stable baseline version (featuring 1D stationary and modal analysis) is located in the **`version0.0`** branch. 
*Please switch to the `version0.0` branch and read its `README.md` for details on features, required libraries, and compilation instructions.*

## Current Development Goals
Active development is taking place to improve the architecture and capabilities of the solver. The current focus includes:
* **Object-Oriented Refactoring:** Restructuring the codebase to be more strictly object-oriented and modular.
* **Optimized Solvers:** Implementing specialized and highly optimized solvers tailored for various classical equations.
* **Boundary Conditions:** Expanding support to handle all remaining types of boundary conditions.

## Future Roadmap
Looking ahead, the project aims to implement:
* **Graphical User Interface (GUI):** To make the solver more accessible and easier to use.
* **Multidimensional Support:** Extending the physics and mathematical models to solve 2D and 3D problems.
=======
**1D FEM Solver - Version 0.1 Update**  
This repository contains a 1D Finite Element Method (FEM) solver. The current update (v0.1) introduces a major architectural overhaul, improving code modularity, usability, and the range of supported physical problems.  
**What's New**  
**1. Object-Oriented Architecture Refactoring**  
The entire codebase has been rewritten to follow strict Object-Oriented Programming (OOP) principles. The code is now much cleaner, highly modular, and easier to maintain or extend in the future.  
**2. CMake Integration**  
Added a CMakeLists.txt file to streamline the build process. The project is now fully compatible with CMake, ensuring an easy, cross-platform compilation process.  
**3. Architecture Separation: **MathSolver ** and **Physics  
To enforce a clear separation of concerns, the core logic has been split into dedicated classes:  
- **MathSolver** ** Class:** A new module exclusively responsible for handling low-level algebraic operations, matrix transformations, and solving linear systems.  
- **Physics** ** Class:** A dedicated class that defines the physical equations governing the problem. This fully decouples the mathematical engine from the physical context.  
**4. User-Friendly API**  
- **Encapsulated Assemblers:** The Matrix_assembler and Vector_assembler methods have been moved to the private scope of the Solver class. This hides the internal mechanics from the end-user, creating a much simpler, robust, and clean public API.  
- **Boundary Conditions Setup:** The definition and assignment of Boundary Conditions have been moved directly to main.cpp. This allows users to set up their simulation parameters quickly without modifying the core library files.  
**5. Advanced Solver Features**  
- **Neumann Boundary Conditions:** The solver now fully supports Neumann (flux/derivative) boundary conditions alongside standard Dirichlet conditions.  
- **Nonlinear Stationary Solver (WIP):** Added support for nonlinear stationary problems. The current implementation utilizes the Newton-Raphson iteration methodto resolve nonlinearities. *Note: This feature is actively being developed (Work in Progress).*  
**6. Quality of Life Improvements**  
- **Automated File Management:** The application now automatically creates the necessary directory tree for storing output results and temporary simulation files. No manual folder  


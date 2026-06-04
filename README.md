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

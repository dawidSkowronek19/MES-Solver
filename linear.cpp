#include "Solver.h"
#include <filesystem>
#include <cstdlib>
#include <string>

int main()
{
    std::cout<<std::string(60, '*')<<"\n\n";
    std::cout<<"\t\t ONE DIMENTION FEM SOLVER\n\n";
    std::cout<<std::string(60, '*')<<"\n";

    std::string compute_mode = "stationary";
    std::string file_dir = "./outdir_"+compute_mode;
    std::filesystem::create_directories(file_dir);

    Grid_1D mesh1(0.0, 1.0, 0.5, 10.0, 12, 200);

    std::vector<BoundaryCondition> bc;
    bc.push_back({0, 1, 1.0});
    bc.push_back({200, 1, 0.5});
    //bc.push_back({150, 2, 0.0});
    mesh1.set_boundaryConditions(bc);
    mesh1.buildGrid();
    mesh1.saveGrid(file_dir);

    MathSolver Math;
    ShapeFunction linear(2);
    //GeneralPDE PDE(Math);
    GeneralSymetricPDE PDE_SYM(Math);
    //Poisson ps(Math);

    Solver solution1(mesh1, linear, PDE_SYM);
    //solution1.Eigen_1D();
    solution1.stationary_1D_implicit();
    solution1.saveSolution(file_dir, compute_mode);
    /*
   Solver EigenSolution(mesh1, linear, PDE);
   EigenSolution.Matrix_assembler("eigen");
   EigenSolution.boundaryConditions("eigen");
   EigenSolution.Eigen_1D();
   EigenSolution.saveSolution("outdir", "eigen");
   */

    int output_py=std::system("python3 ./plot.py");
    return 0;
}
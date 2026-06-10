#include "Solver.h"
#include <filesystem>
#include <cstdlib>
#include <string>

double u(double x) {
    return x<0.5 ? 0.0 : 1.0;}
double v(double x) {return 0.0;}

int main()
{
    std::cout<<std::string(60, '*')<<"\n\n";
    std::cout<<"\t\t ONE DIMENTION FEM SOLVER\n\n";
    std::cout<<std::string(60, '*')<<"\n";

    std::string file_dir;
    std::string compute_mode;

    Grid_1D mesh(0.0, 1.0, 0.5, 10.0, 12, 20);
    std::vector<BoundaryCondition> bc;
    bc.push_back({0, 1, 0.0});
    //bc.push_back({200, 1, 0.0});
    //bc.push_back({150, 2, 0.0});

    mesh.set_boundaryConditions(bc);
    mesh.buildGrid();
    

    MathSolver Math;
    ShapeFunction s_func(4);
    GeneralSymetricPDE PDE_SYM(Math);

    
    //============= EIGEN PART ============
    /*
    compute_mode = "eigen";
    file_dir = "./outdir_"+compute_mode;
    std::filesystem::create_directories(file_dir);

    mesh.saveGrid(file_dir);
    Solver EigenSolution(mesh, s_func, PDE_SYM);
    EigenSolution.Eigen_1D();
    EigenSolution.saveSolution(file_dir, "eigen");
    int output_py_eig=std::system("python3 ./eigen.py");
    */
    //==========================================


    //============ STATIONARY LINEAR PART ======
    /*
    compute_mode = "stationary";
    file_dir = "./outdir_"+compute_mode;
    std::filesystem::create_directories(file_dir);

    mesh.saveGrid(file_dir);
    Solver StationaryLinearSolution(mesh, s_func, PDE_SYM);
    StationaryLinearSolution.stationary_1D_linear();
    StationaryLinearSolution.saveSolution(file_dir, "stationary");
    int output_py_sl=std::system("python3 ./plot.py");
    */
    // ============= Stationary NONLINEAR PART =====
    /*
    compute_mode = "stationary";
    file_dir = "./outdir_"+compute_mode;
    std::filesystem::create_directories(file_dir);

    mesh.saveGrid(file_dir);
    Solver StationaryNonlinearSolution(mesh, s_func, PDE_SYM);
    StationaryNonlinearSolution.stationary_1D_nonlinear();
    StationaryNonlinearSolution.saveSolution(file_dir, "stationary");
    int output_py_snl=std::system("python3 ./plot.py");
    */
    //============= Stationary TIME DEPENDENT ======
    /*
    compute_mode = "time dependent";
    file_dir = "./outdir_"+compute_mode;
    std::filesystem::create_directories(file_dir);

    mesh.saveGrid(file_dir);
    Solver TimeDependentSolution(mesh, s_func, PDE_SYM);
    TimeDependentSolution.timeDependent_1D_linear();
    int output_py_td=std::system("python3 ./gif.py");


    */
   
   
   

    
    return 0;
}
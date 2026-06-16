#include "Solver.h"
#include <filesystem>
#include <cstdlib>
#include <string>

double u(double x) {
    return x<0.5 ? 0.0 : 1.0;}
double v(double x) {return 0.0;}

int main(int argc, char* argv[])
{
    std::cout<<std::string(60, '*')<<"\n\n";
    std::cout<<"\t\t ONE DIMENTION FEM SOLVER\n\n";
    std::cout<<std::string(60, '*')<<"\n";

    std::string input_file="input.inp";
    if (argc == 2) {
        input_file = argv[1];
    } 
    else if (argc > 2) {
        for (int i = 1; i < argc - 1; ++i) {
            if (std::string(argv[i]) == "-in") {
                input_file = argv[i + 1];
                break;
            }
        }
    }

    Parser parser(input_file);
    ConfigParameters config = parser.getParameters();

    std::cout<<config.shapeFunctionDeg;

    Grid_1D mesh(config);
    mesh.set_boundaryConditions();
    mesh.buildGrid();
    MathSolver math;
    ShapeFunction sh_func(config);
    GeneralSymetricPDE PDE(math, config);

    std::string file_dir;
    if (config.work_type == "STATIONARY_LINEAR")
    {
        file_dir = "./outdir_stationary";
        std::filesystem::create_directories(file_dir);

        mesh.saveGrid(file_dir);
        Solver StationaryLinearSolution(mesh, sh_func, PDE, config);
        StationaryLinearSolution.stationary_1D_linear();
        StationaryLinearSolution.saveSolution(file_dir, "stationary");
        int output_py_sl=std::system("python3 ./plot.py");
    }

    else if (config.work_type=="STATIONARY_NONLINEAR")
    {
        file_dir = "./outdir_stationary";
        std::filesystem::create_directories(file_dir);
        

        mesh.saveGrid(file_dir);
        Solver StationaryNonlinearSolution(mesh, sh_func, PDE, config);
        StationaryNonlinearSolution.stationary_1D_nonlinear();
        StationaryNonlinearSolution.saveSolution(file_dir, "stationary");
        int output_py_sl=std::system("python3 ./plot.py");
    }

    else if (config.work_type=="TIME_DEPENDENT")
    {
        file_dir = "outdir_timeDependent";
        std::filesystem::create_directories(file_dir);

        mesh.saveGrid(file_dir);
        Solver TimeDependentSolution(mesh, sh_func, PDE, config);
        TimeDependentSolution.timeDependent_1D_linear();
        int output_py_td=std::system("python3 ./gif.py");
    }

    else if (config.work_type=="EIGEN")
    {
        file_dir = "./outdir_eigen";
        std::filesystem::create_directories(file_dir);

        mesh.saveGrid(file_dir);
        Solver EigenSolution(mesh, sh_func, PDE, config);
        EigenSolution.Eigen_1D();
        EigenSolution.saveSolution(file_dir, "eigen");
        int output_py_eig=std::system("python3 ./eigen.py");
    }



    /*Grid_1D mesh();
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
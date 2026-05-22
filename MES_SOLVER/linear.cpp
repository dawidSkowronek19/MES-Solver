#include "Grid.h"
#include "Solver.h"

int main()
{

    Grid_1D mesh1(-1.0, 1.0, 0.5, 10.0, 12, 200);
    mesh1.set_boundaryConditions();
    mesh1.buildGrid();
    mesh1.saveGrid("./outdir");

    Solver solution1(mesh1);
    solution1.Matrix_assembler("eigen");
    solution1.Vector_assembler();
    solution1.boundaryConditions("eigen");
    //solution1.stationary_1D();
    solution1.Eigen_1D();
    solution1.saveSolution("outdir", "eigen");
    return 0;
}
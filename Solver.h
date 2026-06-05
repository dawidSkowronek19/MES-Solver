#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <complex>
#include <functional>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_splinalg.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "Grid.h"
#include "ShapeFunction.h"
#include "SolverPhysics.h"

class Solver
{

public:
    Solver(Grid_1D &Grid, ShapeFunction &shapefunction, Physics &physics);
    ~Solver();
    void saveSolution(std::string outdir, std::string work_type);
    void stationary_1D_linear();
    void stationary_1D_nonlinear();

    void timeDependent_1D_linear();
    void Eigen_1D();

private:
    Grid_1D &m_grid;
    ShapeFunction &m_shapefunction;
    Physics &m_physics;
    std::vector<double> m_Q, m_U;

    std::vector<double> m_S1D;
    std::vector<double> m_F1D;

    std::vector<double> m_M1D;
    std::vector<double> m_eigenValues1D;
    std::vector<double> m_eigenVectors1D;

    std::vector<double> m_St1D;
    std::vector<double> m_R1D;

    int m_integrationOrder;
    double m_dxSave;
    int m_shapeFunctionDeg;
    std::string m_work_type;

    // integrator
    gsl_integration_fixed_workspace *m_work;
    double *m_local_nodes;
    double *m_local_weights;

    // solver essentials
    static double retOne(double x, double t);
    void local_1D_linear_StifnessMatrix(std::vector<double> &S, int m);
    void local_1D_linear_LoadVector(std::vector<double> &F, int m);
    void local_1D_nonlinear_StifnessMatrix(std::vector<double> &S, int m);
    void local_1D_nonlinear_LoadVector(std::vector<double> &F, int m);
    void local_1D_MassMatrix(std::vector<double> &M, int m, double (*f)(double, double), double t);
    void local_1D_MassMatrix(std::vector<double> &M, int m);
    void local_1D_StiffnessTangentMatrix(std::vector<double> &S_T, int m);

    // assemblers
    void Matrix_assembler(std::string work_type, std::vector<double> &matrix, double (*f)(double x, double t)=nullptr, double t=0.0);
    void Vector_assembler(std::string work_type);

    // boundary conditions
    void boundaryConditions(std::string work_type);

    // solution

    double U_1D(double x);
    double U_1D(double x, int mode);
};

#endif
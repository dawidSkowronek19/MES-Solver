#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <complex>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_splinalg.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "Grid.h"

class Solver{

    public:
        Solver(Grid_1D &grid);
        ~Solver();
        void saveSolution(std::string outdir, std::string work_type);

        void Matrix_assembler(std::string work_type="");
        void Vector_assembler();

        void boundaryConditions(std::string work_type);
        void stationary_1D();
        void Eigen_1D();
    private:
        Grid_1D m_grid;
        std::vector<double> m_Q, m_U;

        std::vector<double> m_S1D;
        std::vector<double> m_F1D;

        std::vector<double> m_M1D;
        std::vector<double> m_eigenValues1D;
        std::vector<double> m_eigenVectors1D;


        int m_integrationOrder;
        double m_dxSave;
        int m_shapeFunctionDeg;

        //integrator
        gsl_integration_fixed_workspace *m_work;
        double *m_local_nodes;
        double *m_local_weights;

        //solver essentials
        void local_1D_StifnessMatrix(std::vector<double> &S, int m);
        void local_1D_LoadVector(std::vector<double> &F, int m);
        void local_1D_MassMatrix(std::vector<double> &M, int m);

        //shape functions
        double Rphi_1D(double ksi, int m);
        double Rphi_1D_deriv(double ksi, int m);

        //algebra solver
        std::vector<double> real_thomas(std::vector<double> &a, std::vector<double> &b, std::vector<double>& c, std::vector<double> &d);


        //solution

        double U_1D(double x);
        double U_1D(double x, int mode);

        double B(double x);
        double C(double x);
        double D(double x);

};

#endif

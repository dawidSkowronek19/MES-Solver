#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_splinalg.h>

struct BoundaryCondition{
    int m;
    int bct; // bounadry condition type 1-Dirichlet, 2-Neumann, 3-Robin
    double bc_value;
};


class Grid_1D{

    public:
        Grid_1D(double x_start, double x_end, double growFactor, double A, int int_order, int numbOfElements);
        void set_boundaryConditions(std::vector<BoundaryCondition>);
        void buildGrid();
        void saveGrid(std::string outdir);

        double get_NodPosition(int m) const;
        double get_elementNumber() const;
        std::vector<BoundaryCondition> get_BC() const;
        std::pair<double, double> get_RodInfo() const;
        ~Grid_1D();

    private:
        std::vector<double> m_nodes;
        std::vector<BoundaryCondition> m_bc;
        double m_x_start, m_x_end;
        double m_growFactor;
        int m_numbOfElements;
        int m_N_int;
        double m_A;

};

#endif
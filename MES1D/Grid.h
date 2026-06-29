#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_splinalg.h>
#include "Parser.h"




class Grid_1D{

    public:
        Grid_1D(ConfigParameters &config);
        void set_boundaryConditions();
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
        ConfigParameters m_config;

};

#endif
#ifndef PARSER_H
#define PARSER_H

    #include <iostream>
    #include <vector>
    #include <string>
    #include <map>
    #include <fstream>
    #include <algorithm>
    #include <cctype>
    #include <sstream>
    #include <functional>
    #include <memory>
    #include "exprtk.hpp"

    struct BoundaryCondition{
    int m;
    int bct; // bounadry condition type 1-Dirichlet, 2-Neumann, 3-Robin
    double bc_value;
};
    
    struct ConfigParameters{
        //GRID
        double x_start;
        double x_end;
        double growFactor;
        double A;
        int integrationOrderGrid;
        int numbOfElements;
        //BOUNDARY CONDITION
        std::vector<BoundaryCondition> boundaryConditions;
        //SOLVER 
        std::string work_type;
        double beta;
        double gamma;
        double epsilon;
        double omega;
        double t_max;
        double dt;
        double du;
        double dxSave;
        //SHAPE FUNCTIONS
        int shapeFunctionDeg;
        //Equations
        std::vector<std::function<double(double, double, double)>> Functions;
        //Initial Values
        std::function<double(double)> u;
        std::function<double(double)> v;

    };

    class Parser{
        public:
            Parser(std::string filename);
            ConfigParameters getParameters();
        private:
            std::map<std::string, std::string> m_rawGridParameters;
            std::map<std::string, std::string> m_rawSolverParameters;
            std::map<std::string, std::string> m_rawBoundaryConditions;
            std::map<std::string, std::string> m_rawFunctions;
            std::map<std::string, std::string> m_rawInitialValues;

            std::string m_filename;
            ConfigParameters param;
            std::vector<BoundaryCondition> m_BoundaryCondition;

            void readRawData();
            void convertToParameters();
            std::string trimBlank(const std::string &str);
    };

#endif
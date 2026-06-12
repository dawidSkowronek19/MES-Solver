#ifndef PARSER_H
#define PARSER_H

    #include <iostream>
    #include <vector>
    #include <string>
    #include <map>
    #include <fstream>
    #include <algorithm>
    #include <cctype>

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

    };

    class Parser{
        public:
            Parser::Parser(std::string filename);
        private:
            std::map<std::string, std::string> m_rawGridParameters;
            std::map<std::string, std::string> m_rawSolverParameters;
            std::map<std::string, std::string> m_rawBoundaryConditions;

            std::string m_filename;

            void readRawData();
            void convertToParameters();
            std::string trimBlank(const std::string &str);
    };

#endif
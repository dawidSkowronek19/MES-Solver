#ifndef GRID
#define GRID

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <functional>
#include <cmath>
#include <set>
#include "Point.hpp"



struct Triangle{
    int p1, p2, p3;
};

struct BC{
    double value;
    int type;
    double t_start, t_end;
    std::function <Position(double)> r;
};

struct GeometryParameters{
    int EdgeNodesNumber;
    double h;
    double toleranceRadius;
    int iterationMax;
    double alpha=0.5;
    
    //std::vector<BC> boundaryConditions;
};

class Grid2D{
    public:
        Grid2D(GeometryParameters geoInit);
        void triangular();
        void saveTrianglesPoints();
        void savePointsList();
        void set_boundaryConditions();
        void create_neighboursList();
        void relaxGrid();
    private:
        GeometryParameters m_geoInit;
        
        double m_edgeLen;
        double m_edgeStart, m_edgeEnd;
        std::function<Position(double)> m_edge;

        std::vector<Position> m_pointsList;
        //std::vector<int> m_boundaryConditionIndexes;
        std::vector<Triangle> m_trianglesList;
        std::vector<std::vector<int>> m_neighboursList;
    };

#endif
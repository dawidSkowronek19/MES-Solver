#ifndef GRID
#define GRID

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <functional>
#include <cmath>

struct Position{
    double x, y;
};


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
    Position r_start, r_end;
    int triangles_nbX;
    int triangles_nbY;

    std::vector<BC> boundaryConditions;
};

class Grid2D{
    public:
        Grid2D(GeometryParameters geoInit);
        void triangular();
        void saveTrianglesPoints();
        void savePointsList();
        void set_boundaryConditions();
    private:
        GeometryParameters m_geoInit;
        double m_Lx, m_Ly;

        std::vector<Position> m_pointsList;
        std::vector<BC> m_boundaryConditionIndexes;
        std::vector<Triangle> m_trianglesList;
    };

#endif
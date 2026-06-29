#ifndef GRID
#define GRID

#include <vector>
#include <iostream>
#include <string>
#include <fstream>

struct Position{
    double x, y;
};

std::ostream& operator<<(std::ostream& os, const Position& pos);

struct Triangle{
    Position P1, P2, P3;
};

struct GeometryParameters{
    Position r_start, r_end;
    int triangles_nbX;
    int triangles_nbY;

};

class Grid2D{
    public:
        Grid2D(GeometryParameters geoInit);
        void triangular();
        void saveTriangles();
    private:
        GeometryParameters m_geoInit;
        double m_Lx, m_Ly;

        
        std::vector<Triangle> m_trianglesList;
    };

#endif
#ifndef DOFHANDLER_HPP
#define DOFHANDLER_HPP

#include "Grid.hpp"
#include "Point.hpp"
struct Node{
    int ID;
    Position r;
};

class DoFHandler{

    public:
        DoFHandler(const Grid2D &grid, const int p);
        void countNodes();
    private:
        std::vector<Node> m_nodes;
        int m_p;
        const Grid2D &m_grid;

};


#endif
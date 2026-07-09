#ifndef DOFHANDLER_HPP
#define DOFHANDLER_HPP

#include "Grid.hpp"
#include "Point.hpp"
#include <map>
#include <algorithm>
class DoFHandler{

    public:
        DoFHandler(const Grid2D &grid, const int p);
        void countNodes();

        const std::vector<std::vector<int>>& get_nodesID() const;
        int get_totalDOF() const;
    private:
        std::vector<std::vector<int>> m_nodesID;
        int m_p;
        const Grid2D &m_grid;
        int m_total_dofs;

};

#endif
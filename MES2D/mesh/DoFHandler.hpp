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
        std::map<int, double> get_boundary_dofs(const std::vector<int> &boundaryNodes, std::function<double(Position)>& bc_fc);
    private:
        std::vector<std::vector<int>> m_nodesID;
        std::map<std::pair<int,int>, int> m_edge_to_dof;
        int m_p;
        const Grid2D &m_grid;
        int m_total_dofs;

};

#endif
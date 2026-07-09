#include "DoFHandler.hpp"

DoFHandler::DoFHandler(const Grid2D &grid, const int p) : m_grid(grid), m_p(p) {}

void DoFHandler::countNodes()
{
    int elementNb = m_grid.getElementNb();
    int pointNb = m_grid.getPointNb();

    m_nodesID.resize(elementNb);

    for (int el_idx=0; el_idx<elementNb; el_idx++)
    {
        Triangle e = m_grid.getTriangle(el_idx);
        m_nodesID[el_idx].push_back(e.p1);
        m_nodesID[el_idx].push_back(e.p2);
        m_nodesID[el_idx].push_back(e.p3);
    }

    m_total_dofs = pointNb;
    if (m_p>1)
    {
        
        for (int el_idx=0; el_idx<elementNb; el_idx++)
        {
            Triangle e = m_grid.getTriangle(el_idx);

            std::vector<std::pair<int,int>> edges = {
                {std::min(e.p1, e.p2), std::max(e.p1, e.p2)},
                {std::min(e.p2, e.p3), std::max(e.p2, e.p3)},
                {std::min(e.p1, e.p3), std::max(e.p1, e.p3)}
            };

            for (const auto& edge : edges)
            {
                if (m_edge_to_dof.find(edge)==m_edge_to_dof.end())
                {
                    m_edge_to_dof[edge] = m_total_dofs;
                    m_total_dofs+=(m_p-1);
                }

                int start_dof = m_edge_to_dof[edge];
                for (int p_idx=0; p_idx<m_p-1; p_idx++)
                {
                    m_nodesID[el_idx].push_back(start_dof+p_idx);
                }
            }
        }
    }

    if (m_p>2)
    {
        int dof_inside_el = (m_p-1)*(m_p-2)/2;

        for (int el_idx=0; el_idx<elementNb; el_idx++)
        {
            for (int in_idx=0; in_idx<dof_inside_el; in_idx++)
                m_nodesID[el_idx].push_back(m_total_dofs++);

        }
    }

}


const std::vector<std::vector<int>>& DoFHandler::get_nodesID() const {return m_nodesID;}
int DoFHandler::get_totalDOF() const {return m_total_dofs;}
void DoFHandler::set_boundary_dofs(const int boundaryNodesNb, const std::function<double(Position)> &bc_fc)
{

    for (int bc_idx=0; bc_idx<boundaryNodesNb; bc_idx++)
    {
        Position r_bc = m_grid.getPoint(bc_idx);
        double value = bc_fc(r_bc);
        m_bc_map[bc_idx] = value;
    }

    for (const auto& [edge, start_dof] : m_edge_to_dof)
    {
        int v1 = edge.first;
        int v2 = edge.second;

        if (v1 < boundaryNodesNb && v2 <boundaryNodesNb)
        {
            Position p1 = m_grid.getPoint(v1);
            Position p2 = m_grid.getPoint(v2);
            Position dr = (p2-p1)/m_p;

            for (int p_idx=0; p_idx<m_p-1; p_idx++)
            {
                int current_dof = start_dof + p_idx;
                Position r_bc = p1 + (p_idx+1)*dr;
                double value = bc_fc(r_bc);

                m_bc_map[current_dof] = value;
            }
        }
    }
}

const std::map<int, double>& DoFHandler::get_boundary_dofs() const {return m_bc_map;}
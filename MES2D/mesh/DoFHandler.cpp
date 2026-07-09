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
        std::map<std::pair<int,int>, int> edge_to_dof;

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
                if (edge_to_dof.find(edge)==edge_to_dof.end())
                {
                    edge_to_dof[edge] = m_total_dofs;
                    m_total_dofs+=(m_p-1);
                }

                int start_dof = edge_to_dof[edge];
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
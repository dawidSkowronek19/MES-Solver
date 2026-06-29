#include "Grid.hpp"

// ========== CONSTRUCTOR ============

Grid2D::Grid2D(GeometryParameters geoInit) : m_geoInit(geoInit)
{

    m_Lx = m_geoInit.r_end.x - m_geoInit.r_start.x;
    m_Ly = m_geoInit.r_end.y - m_geoInit.r_start.y;
    std::cout<<std::string(60, '=');
                
    std::cout<<"\n\n\n \t\t\tGRID SECTION \n\n\n";
    std::cout<<"\tx_start= "<<m_geoInit.r_start.x<<" x_end= "<<m_geoInit.r_end.x<<"\n";
    std::cout<<"\tx_end= "<<m_geoInit.r_start.y<<" y_end= "<<m_geoInit.r_end.y<<"\n";
    std::cout<<"\tLx= "<<m_Lx<<" Ly= "<<m_Ly<<"\n";
    
}

//=========== triangular =============

void Grid2D::triangular()
{
    std::vector<Position> pointsList;
    pointsList.resize(m_geoInit.triangles_nbX*m_geoInit.triangles_nbY);

    double dx = m_Lx/(m_geoInit.triangles_nbX-1);
    double dy = m_Ly/(m_geoInit.triangles_nbY-1);
    for (int y_idx=0; y_idx<m_geoInit.triangles_nbY; y_idx++)
    {
        for (int x_idx=0; x_idx<m_geoInit.triangles_nbX; x_idx++)
        {
            int l = x_idx + y_idx*m_geoInit.triangles_nbX;

            pointsList[l] = {m_geoInit.r_start.x + x_idx*dx,
                                m_geoInit.r_start.y + y_idx*dy
                                };
        }
    }

    for (int y_idx=0; y_idx<m_geoInit.triangles_nbY-1; y_idx++)
    {
        for (int x_idx=0; x_idx<m_geoInit.triangles_nbX-1; x_idx++)
        {
            int l = x_idx + y_idx*m_geoInit.triangles_nbX;

            Position DL = pointsList[l];
            Position DR = pointsList[l+1];
            Position UL = pointsList[l+m_geoInit.triangles_nbX];
            Position UR = pointsList[l+1+m_geoInit.triangles_nbX];

            Triangle t1 = {DL, DR, UR};
            Triangle t2 = {DL, UR, UL};

            m_trianglesList.push_back(t1);
            m_trianglesList.push_back(t2);
        }
    }
    
}

// ========== save grid ==========
std::ostream& operator<<(std::ostream& os, const Position& pos)
{
    os<<"("<<pos.x<<","<<pos.y<<")";
    return os;
}
void Grid2D::saveTriangles()
{

    std::ofstream triangle_file("./outdir/trianglesList.dat");

    for (const auto& triangle : m_trianglesList)
    {
        triangle_file << triangle.P1<<" "<<triangle.P2<<" "<<triangle.P3<<"\n";
    }


}
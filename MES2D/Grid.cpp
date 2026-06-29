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
    m_pointsList.resize(m_geoInit.triangles_nbX*m_geoInit.triangles_nbY);

    double dx = m_Lx/(m_geoInit.triangles_nbX-1);
    double dy = m_Ly/(m_geoInit.triangles_nbY-1);
    for (int y_idx=0; y_idx<m_geoInit.triangles_nbY; y_idx++)
    {
        for (int x_idx=0; x_idx<m_geoInit.triangles_nbX; x_idx++)
        {
            int l = x_idx + y_idx*m_geoInit.triangles_nbX;

            m_pointsList[l] = {m_geoInit.r_start.x + x_idx*dx,
                                m_geoInit.r_start.y + y_idx*dy
                                };
        }
    }

    for (int y_idx=0; y_idx<m_geoInit.triangles_nbY-1; y_idx++)
    {
        for (int x_idx=0; x_idx<m_geoInit.triangles_nbX-1; x_idx++)
        {
            int l = x_idx + y_idx*m_geoInit.triangles_nbX;

            int DL = l;
            int DR = l+1;
            int UL = l+m_geoInit.triangles_nbX;
            int UR = l+1+m_geoInit.triangles_nbX;

            

            Triangle t1 = {DL, DR, UR};
            Triangle t2 = {DL, UR, UL};

            m_trianglesList.push_back(t1);
            m_trianglesList.push_back(t2);
        }
    }
    
}

void Grid2D::set_boundaryConditions()
{
    //for test 
    // ********************************************
    int numb_of_points_onBC=10;
    m_pointsList.resize(numb_of_points_onBC);
    double radius=0.2;
    double x_0 = (m_geoInit.r_end.x - m_geoInit.r_start.x)/2.0;
    double y_0 = (m_geoInit.r_end.y - m_geoInit.r_start.y)/2.0;

    auto r_func =[radius, x_0, y_0](double t)->Position{
        return {x_0+radius*cos(t), y_0+radius*sin(t)};
    };
    m_geoInit.boundaryConditions.push_back({0.0, 1,0.0, 2.0*M_PI, r_func});
    // ********************************************


    for (const auto &boundaryCondition : m_geoInit.boundaryConditions)
    {
        double dt = (boundaryCondition.t_end - boundaryCondition.t_start)/(numb_of_points_onBC-1);
        for (int idx=0; idx<numb_of_points_onBC; idx++)
        {
            double acc_t = boundaryCondition.t_start + dt*idx;
            Position acc_R = boundaryCondition.r(acc_t);
            m_pointsList[idx] = acc_R;
        }
    }


}

// ========== save grid ==========
void Grid2D::saveTrianglesPoints()
{

    std::ofstream triangle_file("./outdir/trianglesPointsList.dat");

    for (const auto& triangle : m_trianglesList)
    {
        triangle_file << triangle.p1<<" "<<triangle.p2<<" "<<triangle.p3<<"\n";
    }


}

void Grid2D::savePointsList()
{
    std::ofstream point_lists_file("./outdir/pointsList.dat");

    for (const auto& point : m_pointsList)
    {
        point_lists_file<<point.x<<" "<<point.y<<"\n";
    }

}
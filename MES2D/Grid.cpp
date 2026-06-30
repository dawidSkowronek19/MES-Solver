#include "Grid.hpp"

// ========== CONSTRUCTOR ============

Grid2D::Grid2D(GeometryParameters geoInit) : m_geoInit(geoInit)
{

    std::cout<<std::string(60, '=');
    //for tests

    auto edge_func = [](double t)->Position{
        return {cos(t),sin(t)};
    };

    m_edge = edge_func;
    m_edgeStart=0.0;
    m_edgeEnd=2.0*M_PI;

    std::cout<<"\n\n\n \t\t\tGRID SECTION \n\n\n";
    //std::cout<<"\tL= "<<m_edgeLen<<"\n";
    
}

//=========== triangular =============

void Grid2D::triangular()
{
    std::cout<<"EDGE POINTS... "<<std::flush;
    double dt = (m_edgeEnd-m_edgeStart)/(m_geoInit.EdgeNodesNumber-1);
    std::vector<Position> oldRing;
    std::vector<Position> newRing;

    for (int point_idx=0; point_idx<m_geoInit.EdgeNodesNumber; point_idx++)
    {
        double t = m_edgeStart + dt*point_idx;
        Position point = m_edge(t);
        oldRing.push_back(point);
        m_pointsList.push_back(point);
    }
    std::cout<<"DONE\n";

    auto calculate_normals = [&oldRing](int idx)->Position{
        Position A = oldRing[idx-1];
        Position C = oldRing[idx+1];

        Position CA = {C.x-A.x, C.y-A.y};

        return {-CA.y, CA.x};
    };
    std::vector<Position> newRing;

    while(true)
    {
        Position n = calculate_normals(0);
        Position newPoint = {oldRing[0].x + n.x, oldRing[0].y + n.y};
        newRing.push_back(newPoint);

        for (int idx=1; idx<oldRing.size(); idx++)
        {
            n = calculate_normals(idx);
            newPoint = {oldRing[idx].x + n.x, oldRing[idx].y + n.y};
            
            double dist2 = (newPoint.x - newRing[idx-1].x)*(newPoint.x - newRing[idx-1].x) + (newPoint.y - newRing[idx-1].y)*(newPoint.y - newRing[idx-1].y);
            
            if (dist2<m_geoInit.toleranceRadius)
            {
                Position avg;
                Position previous = newRing[idx-1];
                
                avg.x = (newPoint.x + previous.x)/2.0;
                avg.y = (newPoint.y +previous.y)/2.0;

                newRing.pop_back();
                newRing.push_back(avg);
            }
            else
                newRing.push_back(newPoint);

            if (newRing.size()==0)
                break;
        }

        
        for (int idx=0; idx<oldRing.size()-1; idx++)
        {
            Position BL = oldRing[idx];
            Position BR = oldRing[idx+1];
        }

        oldRing=newRing;
        newRing.clear();

    }
    
}

/*void Grid2D::set_boundaryConditions()
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
            m_boundaryConditionIndexes.push_back(idx);
        }
    }


}*/

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
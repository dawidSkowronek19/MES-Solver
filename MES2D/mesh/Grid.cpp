#include "Grid.hpp"

// ========== CONSTRUCTOR ============

Grid2D::Grid2D(GeometryParameters geoInit) : m_geoInit(geoInit)
{

    std::cout<<std::string(60, '=');
    //for tests

    auto edge_func = [](double t)->Position{
        double x, y;
        /*if (t<=0.25)
        {
            x=t;
            y=0.0;
        }
        else if (t>0.25 && t<=0.5)
        {
            x=0.25;
            y=t-0.25;
        }
        else if (t>0.5 && t<=0.75)
        {
            x=0.75-t;
            y=0.25;
        }
        else
        {
            x=0.0;
            y=1.0-t;
        }
        */
        
        x=0.7*cos(2*M_PI*t);
        y=sin(2*M_PI*t);


        return {x,y};
    };
    

    m_edge = edge_func;
    m_edgeStart=0.0;
    m_edgeEnd=1.0;

    std::cout<<"\n\n\n \t\t\tGRID SECTION \n\n\n";
    //std::cout<<"\tL= "<<m_edgeLen<<"\n";
    
}

//=========== triangular =============

void Grid2D::triangular()
{
    std::cout<<"EDGE POINTS... "<<std::flush;
    double dt = (m_edgeEnd-m_edgeStart)/(m_geoInit.EdgeNodesNumber);
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

    auto calculate_normals = [&oldRing, this](int idx)->Position{
        int size = oldRing.size();
        Position A = oldRing[(idx - 1 + size) % size];
        Position B = oldRing[idx];
        Position C = oldRing[(idx + 1) % size];

        Position BA = B-A;
        Position normal1 = {-BA.y(), BA.x()};
        normal1/=normal1.len();

        Position CB = C-B;
        Position normal2 = {-CB.y(), CB.x()};
        normal2/=normal2.len();


        Position normalB = normal1 + normal2;
        normalB/=normalB.len();


        double dot = normal1.dot(normal2);
        double scale = m_geoInit.h*sqrt(2.0/(1.0+dot));
    

        return normalB*scale;
    };

    while(true)
    {
        
        
        static int iter = 0;
        std::cout<<oldRing.size()<<"\n";
        if (oldRing.size()<3 || ++iter>10000)
            break;

        Position avg={0.0, 0.0};
        for (const auto& p : oldRing) { 
            avg+=p; 
        }

        avg /=oldRing.size();

        double max_dist2 = 0.0;
        for (const auto& p : oldRing) 
        {
            Position diff = p-avg;
            double dist2 = diff.len2();

            if (dist2 > max_dist2) max_dist2 = dist2;
        }


        if (max_dist2 <= m_geoInit.h * m_geoInit.h * 1.5)
        {

            int oldRingStartIdx = m_pointsList.size() - oldRing.size();
            

            m_pointsList.push_back(avg);
            int centerIdx = m_pointsList.size() - 1;


            for (int idx = 0; idx < oldRing.size(); idx++) 
            {
                int globBL = oldRingStartIdx + idx;
                int globBR = oldRingStartIdx + ((idx + 1) % oldRing.size());
                m_trianglesList.push_back({globBL, globBR, centerIdx});
            }
            
            std::cout << "MESH CLOSED\n";
            break;
        }
        // ==========================================

        std::vector<int> mapping(oldRing.size());


        Position n = calculate_normals(0);
        Position newPoint = oldRing[0] + n;
        newRing.push_back(newPoint);
        mapping[0]=0;

        for (int idx=1; idx<oldRing.size(); idx++)
        {
            n = calculate_normals(idx);
            newPoint = oldRing[idx] + n;

            //sewwing logic
            double dist2 = (newPoint - newRing.back()).len2();
            if (dist2>=m_geoInit.toleranceRadius*m_geoInit.toleranceRadius)
                newRing.push_back(newPoint);
            mapping[idx] = newRing.size()-1;
        }

        double dist2LastFirst = (newRing.front()-newRing.back()).len2();
        if (newRing.size()>1 && dist2LastFirst<m_geoInit.toleranceRadius*m_geoInit.toleranceRadius)
        {   
            newRing.pop_back();
            for (int i=0; i<mapping.size(); i++)
            {
                if (mapping[i] == newRing.size())
                    mapping[i]=0;
            }
        }
        int newRing_PointsNumber = newRing.size();
        int oldRing_PointsNumber = oldRing.size();
        int pointList_PointsNumber = m_pointsList.size();
        m_pointsList.reserve(pointList_PointsNumber+newRing_PointsNumber);

        for (int idx=0; idx<newRing_PointsNumber; idx++)
        {
            m_pointsList.push_back(newRing[idx]);
        }


        //creating triangles      
        int previous_step_pointList_PointsNumber = pointList_PointsNumber - oldRing_PointsNumber;  
        for (int idx=0; idx<oldRing.size(); idx++)
        {
            int loc_BL = idx; //old ring
            int loc_BR = (idx + 1) % oldRing.size(); //old ring
            int loc_TL = mapping[loc_BL]; //new ring
            int loc_TR = mapping[loc_BR]; //new ring


            int globBL = previous_step_pointList_PointsNumber +loc_BL;
            int globBR = previous_step_pointList_PointsNumber +loc_BR;

            int globTL = pointList_PointsNumber + loc_TL;
            int globTR = pointList_PointsNumber + loc_TR;
            if (loc_TR!=loc_TL)
            {
                if (idx % 2 == 0) 
                {
                    m_trianglesList.push_back({globBL, globBR, globTR});
                    m_trianglesList.push_back({globBL, globTR, globTL});
                } 
                else 
                {
                    m_trianglesList.push_back({globBL, globBR, globTL});
                    m_trianglesList.push_back({globBR, globTR, globTL});
                }
            }
            else //for degeneration
            {
                m_trianglesList.push_back({globBL, globBR, globTR});
            }

        }

        oldRing=newRing;
        newRing.clear();
        newRing.reserve(oldRing_PointsNumber);

    }
    
}

void Grid2D::create_neighboursList()
{
    std::vector<std::set<int>> neighboursList(m_pointsList.size());
    m_neighboursList.resize(m_pointsList.size());

    for (const auto triangle:m_trianglesList)
    {
        neighboursList[triangle.p1].insert(triangle.p2);
        neighboursList[triangle.p1].insert(triangle.p3);

        neighboursList[triangle.p2].insert(triangle.p1);
        neighboursList[triangle.p2].insert(triangle.p3);

        neighboursList[triangle.p3].insert(triangle.p1);
        neighboursList[triangle.p3].insert(triangle.p2);

    }

    for (int idx=0; idx<neighboursList.size(); idx++)
    {
        m_neighboursList[idx].assign(neighboursList[idx].begin(), neighboursList[idx].end());
    }


}

void Grid2D::relaxGrid()
{
    std::vector<Position> buffor=m_pointsList;
    for (int iteration=0; iteration<m_geoInit.iterationMax; iteration++)
    {
        
        for(int idx = m_geoInit.EdgeNodesNumber; idx<m_pointsList.size(); idx++)
        {
            int neighboursNb = m_neighboursList[idx].size();
            Position newPoint = {0.0, 0.0};

            for (const auto& neighbour_idx : m_neighboursList[idx])
            {
                Position neighbour = m_pointsList[neighbour_idx];
                newPoint+=neighbour;
            }

            newPoint = m_geoInit.alpha*newPoint/neighboursNb + (1.0-m_geoInit.alpha)*m_pointsList[idx];

            buffor[idx]=newPoint;
        }

        m_pointsList=buffor;

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
        point_lists_file<<point<<"\n";
    }

}
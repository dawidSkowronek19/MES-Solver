#include "Grid.hpp"



int main()
{
    GeometryParameters geoparam{
        50,
        2*sin(M_PI/50),
        2*0.7*sin(M_PI/50),
        100
    };

    Grid2D testgrid(geoparam);
    testgrid.triangular();
    testgrid.create_neighboursList();
    testgrid.relaxGrid();
    testgrid.savePointsList();
    testgrid.saveTrianglesPoints();

}

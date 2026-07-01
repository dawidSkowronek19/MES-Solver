#include "Grid.hpp"



int main()
{
    GeometryParameters geoparam{
        18,
        2*sin(M_PI/18),
        2*0.3*sin(M_PI/18),
        100
    };

    Grid2D testgrid(geoparam);
    testgrid.triangular();
    testgrid.create_neighboursList();
    testgrid.relaxGrid();
    testgrid.savePointsList();
    testgrid.saveTrianglesPoints();

}

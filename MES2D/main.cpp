#include "./mesh/Grid.hpp"



int main()
{
    GeometryParameters geoparam{
        100,
        3*sin(M_PI/100),
        3*0.8*sin(M_PI/100),
        100
    };

    Grid2D testgrid(geoparam);
    testgrid.triangular();
    testgrid.create_neighboursList();
    testgrid.relaxGrid();
    testgrid.savePointsList();
    testgrid.saveTrianglesPoints();

}

#include "Grid.hpp"



int main()
{
    GeometryParameters geoparam{
        40,
        0.02,
        0.01
    };

    Grid2D testgrid(geoparam);
    testgrid.triangular();
    testgrid.savePointsList();
    testgrid.saveTrianglesPoints();

}

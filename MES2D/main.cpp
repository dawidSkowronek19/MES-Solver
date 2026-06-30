#include "Grid.hpp"



int main()
{
    GeometryParameters geoparam{
        10
    };

    Grid2D testgrid(geoparam);
    testgrid.triangular();
    testgrid.savePointsList();
    testgrid.saveTrianglesPoints();

}

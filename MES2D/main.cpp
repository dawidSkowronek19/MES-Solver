#include "Grid.hpp"



int main()
{
    GeometryParameters geoparam{
        {0.0, 0.0}, {1.0, 2.0},
        10,
        20
    };

    Grid2D testgrid(geoparam);
    testgrid.triangular();
    testgrid.saveTriangles();

}

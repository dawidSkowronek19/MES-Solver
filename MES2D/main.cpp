#include "./mesh/Grid.hpp"
#include "./solver/shapeFunction.hpp"
#include "./solver/Quadrature.hpp"
#include "./solver/Physics.hpp"
#include "./solver/Assembler.hpp"

#include <fstream>

int main()
{
    GeometryParameters geoparam{
        100,
        2*sin(M_PI/100),
        2*0.5*sin(M_PI/100),
        100
    };

}   

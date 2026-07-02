#include "shapeFunction.hpp"

Position ShapeFunction::bar_to_cartes(double ksi, double eta, int idx) const
{
    Triangle triangle = m_grid.getTriangle(idx);
    Position p1 = m_grid.getPoint(triangle.p1);
    Position p2 = m_grid.getPoint(triangle.p2);
    Position p3 = m_grid.getPoint(triangle.p3);
    
    double x = 
}
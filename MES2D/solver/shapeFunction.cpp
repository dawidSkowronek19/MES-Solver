#include "shapeFunction.hpp"

void ShapeFunction::find_accElementPoints(Triangle Element)
{
    m_accElement={
        m_grid.getPoint(Element.p1),
        m_grid.getPoint(Element.p2),
        m_grid.getPoint(Element.p3)
    };
}

Position ShapeFunction::XY(const Position &p1, const Position &p2, const Position &p3) const
{
    double x = p1.x() +(p2.x()-p1.x())*m_ksi + (p3.x() - p1.x())*m_eta;
    double y = p1.y() + (p2.y()-p1.y())*m_ksi + (p3.y() - p1.y())*m_eta;

    return {x,y};
}

Position ShapeFunction::loc_to_cartes(double ksi, double eta)
{
    Position p1 = m_accElement.p1;
    Position p2 = m_accElement.p2;
    Position p3 = m_accElement.p3;

    m_ksi=ksi; m_eta = eta;
    Position r = XY(p1, p2, p3);

    return r;
}

Position ShapeFunction::loc_to_cartes(double ksi, double eta, Triangle element)
{
    Position p1 = m_grid.getPoint(element.p1);
    Position p2 = m_grid.getPoint(element.p2);
    Position p3 = m_grid.getPoint(element.p3);

    m_ksi = ksi; m_eta = eta;
    Position r = XY(p1,p2,p3);

    return r;
}
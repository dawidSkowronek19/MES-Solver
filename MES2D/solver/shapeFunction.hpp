#ifndef SHAPE_HPP
#define SHAPE_HPP
#include "../mesh/Grid.hpp"

struct ElementPointPositions{
    Position p1, p2, p3;
};

class ShapeFunction{

    public:
        double phi(double ksi ,double eta, Triangle element) const;
        std::tuple<Position, Position, Position> element_idx;
        void find_accElementPoints(Triangle Element);
        Position loc_to_cartes(double ksi ,double eta); //with preconditioning, requires find_accElementPoints() before
        Position loc_to_cartes(double ksi, double eta, Triangle element); //

    private:
        const Grid2D& m_grid;
        ElementPointPositions m_accElement;
        double m_ksi, m_eta;

        Position XY(const Position &p1, const Position &p2, const Position &p3) const;


};

#endif
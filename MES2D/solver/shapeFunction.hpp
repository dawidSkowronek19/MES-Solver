#ifndef SHAPE_HPP
#define SHAPE_HPP
#include "../mesh/Grid.hpp"

struct Element{
    int idx;
    Triangle triangle;
    Position p1
};

class ShapeFunction{

    public:
        double phi(double ksi ,double eta, int i) const;
        std::tuple<Position, Position, Position> element_idx;
        Position bar_to_cartes(double ksi, double eta, int idx) const;

    private:
        const Grid2D& m_grid;

};

#endif
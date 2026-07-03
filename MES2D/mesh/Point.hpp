#ifndef POINT_HPP
#define POINT_HPP
#include <iostream>
#include <cmath>
class Position{
    public:
        Position(double x, double y);
        Position();

        double dot(const Position vec) const;
        double len2() const;
        double len() const;
        Position operator*(const double alpha) const;
        Position operator+(const Position &other) const;
        Position operator-(const Position &other) const;
        Position operator/(const double alpha) const;

        Position& operator/=(const double alpha);
        Position& operator+=(const double alpha);
        Position& operator+=(const Position &other);
        Position& operator*=(const double alpha);
        Position& operator-=(const double alpha);

        friend std::ostream& operator<<(std::ostream& os, const Position& pos);
        friend Position operator*(double alpha, const Position& pos);

        //getters
        double x() const;
        double y() const;

    private:
        double m_x;
        double m_y;
};

#endif
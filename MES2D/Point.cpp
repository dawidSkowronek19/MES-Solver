#include "Point.hpp"

std::ostream& operator<<(std::ostream& os, const Position& pos) {
    os << pos.m_x << " " << pos.m_y;
    return os;
}
Position operator*(double alpha, const Position& pos) {
    return pos * alpha; 
}

//================ POSITION CLASS =====================

Position::Position(double x, double y) : m_x(x), m_y(y) {}


double Position::x() const {return m_x;}
double Position::y() const {return m_y;}

double Position::dot(const Position vec) const{
    return m_x*vec.m_x + m_y*vec.m_y;
}

double Position::len2() const{
    return m_x*m_x+m_y*m_y;
}
double Position::len() const{
    return sqrt(len2());
}

Position Position::operator*(const double alpha) const{
    return {alpha*m_x, alpha*m_y};
}

Position Position::operator+(const Position &other) const{
    return {m_x+other.m_x, m_y+other.m_y};
}

Position Position::operator-(const Position &other) const{
    return {m_x-other.m_x, m_y-other.m_y};
}

Position Position::operator/(const double alhpa) const{
    return {m_x/alhpa, m_y/alhpa};
}

Position& Position::operator/=(const double alpha)
{
    m_x/=alpha;
    m_y/=alpha;

    return *this;
}
Position& Position::operator+=(const double alpha)
{
    m_x+=alpha;
    m_y+=alpha;

    return *this;
}
Position& Position::operator+=(const Position &other)
{
    m_x+=other.m_x;
    m_y+=other.m_y;

    return *this;
}
Position& Position::operator*=(const double alpha)
{
    m_x*=alpha;
    m_y*=alpha;

    return *this;
}
Position& Position::operator-=(const double alpha)
{
    m_x-=alpha;
    m_y-=alpha;

    return *this;
}
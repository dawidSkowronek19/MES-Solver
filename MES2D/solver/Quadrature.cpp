#include "Quadrature.hpp"


Quadrature::Quadrature() {}
double Quadrature::get_ksi(const int idx) const {return m_integrationPoints[idx].first;}
double Quadrature::get_eta(const int idx) const {return m_integrationPoints[idx].second;}
double Quadrature::get_weight(const int idx) const {return m_weights[idx];}
int Quadrature::get_point_numb() const {return m_weights.size();}
const std::vector<std::pair<double, double>>& Quadrature::get_integrationPoints() const {return m_integrationPoints;}


void Quadrature::add_point(std::pair<double, double> point, double weight)
{
    m_integrationPoints.push_back(point);
    m_weights.push_back(weight);
}

Quadrature QuadratureFactory::refTriangle(int p)
{
    Quadrature quad;

    if (p == 1) 
    {
        quad.add_point({1.0/3.0, 1.0/3.0}, 1.0/2.0);
    }
    else if (p == 2) 
    {

        quad.add_point({0.091576213, 0.091576213}, 0.109951744);
        quad.add_point({0.091576213, 0.816847574}, 0.109951744);
        quad.add_point({0.816847574, 0.091576213}, 0.109951744);
        quad.add_point({0.445948491, 0.445948491}, 0.223381589);
        quad.add_point({0.445948491, 0.108103018}, 0.223381589);
        quad.add_point({0.108103018, 0.445948491}, 0.223381589);
    }
    else if (p == 3) 
    {
        quad.add_point({1.0/3.0, 1.0/3.0}, -27.0/96.0);
        quad.add_point({1.0/5.0, 1.0/5.0},  25.0/96.0);
        quad.add_point({3.0/5.0, 1.0/5.0},  25.0/96.0);
        quad.add_point({1.0/5.0, 3.0/5.0},  25.0/96.0);
        }
    
    return quad;
}
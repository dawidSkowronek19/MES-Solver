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
        double w = 1.0 / 6.0;
        quad.add_point({1.0/6.0, 1.0/6.0}, w);
        quad.add_point({2.0/3.0, 1.0/6.0}, w);
        quad.add_point({1.0/6.0, 2.0/3.0}, w);
    }
    else if (p == 2) 
    {

        double w1 = 0.109951743655322 / 2.0; 
        double w2 = 0.223381589678011 / 2.0; 

        quad.add_point({0.091576213509771, 0.091576213509771}, w1);
        quad.add_point({0.091576213509771, 0.816847572980459}, w1);
        quad.add_point({0.816847572980459, 0.091576213509771}, w1);
        quad.add_point({0.445948490915965, 0.445948490915965}, w2);
        quad.add_point({0.445948490915965, 0.108103018168070}, w2);
        quad.add_point({0.108103018168070, 0.445948490915965}, w2);
    }
    else if (p == 3) 
    {
        double w1 = 0.025422453185103;
        double w2 = 0.058393137863189;
        double w3 = 0.041425537809186;
        double w4 = 0.041425537809186;

        quad.add_point({0.063089014491502, 0.063089014491502}, w1);
        quad.add_point({0.063089014491502, 0.873821971016996}, w1);
        quad.add_point({0.873821971016996, 0.063089014491502}, w1);
        
        quad.add_point({0.249286745170910, 0.249286745170910}, w2);
        quad.add_point({0.249286745170910, 0.501426509658180}, w2);
        quad.add_point({0.501426509658180, 0.249286745170910}, w2);
        
        quad.add_point({0.310352451033784, 0.053145049844816}, w3);
        quad.add_point({0.053145049844816, 0.310352451033784}, w3);
        quad.add_point({0.310352451033784, 0.636502499121400}, w3);
        
        quad.add_point({0.636502499121400, 0.310352451033784}, w4);
        quad.add_point({0.053145049844816, 0.636502499121400}, w4);
        quad.add_point({0.636502499121400, 0.053145049844816}, w4);
    }
    
    return quad;
}
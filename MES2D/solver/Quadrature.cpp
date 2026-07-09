#include "Quadrature.hpp"


Quadrature::Quadrature() {}
double Quadrature::get_ksi(const int idx) const {return m_integrationPoints[idx].first;}
double Quadrature::get_eta(const int idx) const {return m_integrationPoints[idx].second;}
double Quadrature::get_weight(const int idx) const {return m_weights[idx];}
int Quadrature::get_point_numb() const {return m_weights.size();}
const std::vector<std::pair<double, double>>& Quadrature::get_integrationPoints() const {return m_integrationPoints;}

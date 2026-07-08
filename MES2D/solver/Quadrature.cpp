#include "Quadrature.hpp"


Quadrature::Quadrature() {}
double Quadrature::get_ksi(const int idx) const {return m_ksi[idx];}
double Quadrature::get_eta(const int idx) const {return m_eta[idx];}
double Quadrature::get_weight(const int idx) const {return m_weights[idx];}
int Quadrature::get_point_numb() const {return m_weights.size();}

#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP


#include <vector>

class Quadrature{
    public:
        Quadrature();
        double get_ksi(const int idx) const;
        double get_eta(const int idx) const;
        double get_weight(const int idx) const;
        int get_point_numb() const;

    private:
        std::vector<double> m_ksi;
        std::vector<double> m_eta;
        std::vector<double> m_weights;
};

#endif
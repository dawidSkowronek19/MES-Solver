#ifndef SHAPE_HPP
#define SHAPE_HPP
#include "../mesh/Grid.hpp"
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <tuple>

struct ElementPointPositions{
    Position p1, p2, p3;
};

class ShapeFunction{

    public:
        ShapeFunction(const ElementPointPositions &acc_element, const int p);
        double phi(double ksi ,double eta, int idx) const;
        std::tuple<Eigen::Matrix2d, Eigen::Matrix2d, double> JacobiEssentials() const;
        std::pair<double, double> div_phi(double ksi, double eta, int idx) const;
        Position loc_to_cartes(double ksi ,double eta) const;


    private:
        const int m_p;
        ElementPointPositions m_accElement;
        std::vector<std::tuple<int, int, int>> m_IJK;

        std::tuple<double, double, double> lambda_gen(const double ksi, const double eta) const;
        double Silvester(const int i, const double lambda) const;
        double divSilvester(const int i, const double lambda) const;

        Position KsiEta_to_XY(double ksi, double eta) const;


};

#endif
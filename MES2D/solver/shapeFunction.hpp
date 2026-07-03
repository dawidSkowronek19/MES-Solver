#ifndef SHAPE_HPP
#define SHAPE_HPP
#include "../mesh/Grid.hpp"
#include <vector>
#include <Eigen/Dense>
#include <tuple>

struct ElementPointPositions{
    Position p1, p2, p3;
};

class ShapeFunction{

    public:
        ShapeFunction(Grid2D& grid, int p);
        double phi(double ksi ,double eta, int idx);
        std::tuple<Eigen::Matrix2d, Eigen::Matrix2d, double> JacobiEssentials();
        std::pair<double, double> div_phi(double ksi, double eta, int idx);
        void find_accElementPoints(Triangle Element);
        Position loc_to_cartes(double ksi ,double eta); //with preconditioning, requires find_accElementPoints() before
        Position loc_to_cartes(double ksi, double eta, Triangle element); //

    private:
        const Grid2D& m_grid;
        const int m_p;
        ElementPointPositions m_accElement;
        double m_ksi, m_eta;
        std::vector<std::tuple<int, int, int>> m_IJK;

        std::tuple<double, double, double> lambda_gen();
        double Silvester(const int i, const double lambda) const;
        double divSilvester(const int i, const double lambda) const;

        Position XY(const Position &p1, const Position &p2, const Position &p3) const;


};

#endif
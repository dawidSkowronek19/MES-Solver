#ifndef SHAPE_HPP
#define SHAPE_HPP
#include "../mesh/Grid.hpp"
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <tuple>



struct ShapeFunctionEssentials{
    double phi, dphi_dksi, dphi_deta;
};

class ShapeFunction{

    public:
        ShapeFunction(const int p);
        
        double get_phi(const int idx, const int sh_nb) const;
        double get_dphi_dksi(const int idx, const int sh_nb) const;
        double get_dphi_deta(const int idx, const int sh_nb) const;
        void set_cached_values(const std::vector<std::pair<double, double>>&integration_points);
        int get_p() const;


    private:
        const int m_p;
        std::vector<std::tuple<int, int, int>> m_IJK;
        
        Eigen::MatrixXd m_phi, m_dphi_dksi, m_dphi_deta;

        std::tuple<double, double, double> lambda_gen(const double ksi, const double eta) const;
        double Silvester(const int i, const double lambda) const;
        double divSilvester(const int i, const double lambda) const;
        double phi(double ksi ,double eta, int idx) const;
        std::pair<double, double> div_phi(double ksi, double eta, int idx) const;
        
        //Position KsiEta_to_XY(double ksi, double eta) const;


};


#endif
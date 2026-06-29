#ifndef MATH_H
#define MATH_H

#include <Eigen/Dense>
#include <iostream>

class MathSolver{
    public:
        std::pair<Eigen::VectorXd, Eigen::MatrixXd>EigenSymetric_real(std::vector<double> &S, std::vector<double> &M);
        std::pair<Eigen::VectorXd, Eigen::MatrixXd> EigenGeneral_real(std::vector<double> &S, std::vector<double> &M);

        Eigen::VectorXd StEqSymetric_real(std::vector<double> &S, std::vector<double> &F);
        Eigen::VectorXd StEqGeneral_real(std::vector<double> &S, std::vector<double> &F);


};

#endif
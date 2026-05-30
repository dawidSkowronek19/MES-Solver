#ifndef PHYSICS_H
#define PHYSICS_H

#include <iostream>
#include "MathSolver.h"
class Physics{
    public:
        Physics(MathSolver &math);
        virtual ~Physics() = default;

        
        virtual Eigen::VectorXd solver(std::vector<double> &S, std::vector<double> &F)=0;
        virtual std::pair<Eigen::VectorXd, Eigen::MatrixXd> EigenSolver(std::vector<double> &S, std::vector<double> &M);

        virtual double B(double x, double u);
        virtual double C(double x, double u);
        virtual double D(double x, double u);


    protected:
        MathSolver &m_math;
};

class Poisson : public Physics{

    // d2/dx2 u(x) = D(x)

    public:
        Poisson(MathSolver &math);
        Eigen::VectorXd solver(std::vector<double> &S, std::vector<double> &F) override;
        double D(double x, double u) override;
};


class GeneralPDE : public Physics{
        // [d2/dx2 +B(x) d/dx +C(x)] u(x) = D(x)

        // [d2/dx2 + B(x) d/dx +C(x)]u(x)=\lambda u(x)

    public:
        GeneralPDE(MathSolver &math);
        Eigen::VectorXd solver(std::vector<double> &S, std::vector<double> &F) override;
        std::pair<Eigen::VectorXd, Eigen::MatrixXd> EigenSolver(std::vector<double> &S, std::vector<double> &M) override;

        double B(double x, double u) override;
        double C(double x, double u) override;
        double D(double x, double u) override;

};

class GeneralSymetricPDE : public Physics{
    // [d2/dx2 +C(x)]u(x)=D(x)
    // [d2/dx2 +C(x)]u(x)=\lambda u(x)
    public:
        GeneralSymetricPDE(MathSolver &math);
        Eigen::VectorXd solver(std::vector<double> &S, std::vector<double> &F) override;
        std::pair<Eigen::VectorXd, Eigen::MatrixXd> EigenSolver(std::vector<double> &S, std::vector<double> &M) override;

        double C(double x, double u) override;
        double D(double x, double u) override;

};

#endif
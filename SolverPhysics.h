#ifndef PHYSICS_H
#define PHYSICS_H

#include <iostream>
#include "MathSolver.h"
class Physics{
    public:
        Physics(MathSolver &math);
        virtual ~Physics() = default;

        virtual void connectTime(const double *t);
        
        virtual Eigen::VectorXd solver(std::vector<double> &S, std::vector<double> &F)=0;
        virtual std::pair<Eigen::VectorXd, Eigen::MatrixXd> EigenSolver(std::vector<double> &S, std::vector<double> &M);

        virtual double B(double x, double u);
        virtual double C(double x, double u);
        virtual double D(double x, double u);
        //virtual double E(double x, double u);
        //virtual double F(double x, double u);

        virtual double dB_du(double x, double u);
        virtual double dC_du(double x, double u);
        virtual double dD_du(double x, double u);

        virtual double B(double x);
        virtual double C(double x);
        virtual double D(double x);
        virtual double E(double x);
        virtual double F(double x);
        virtual double getCurrentTime();

    protected:
        MathSolver &m_math;
    private:
        const double *m_t = nullptr;
};

class Poisson : public Physics{

    // d2/dx2 u(x) = D(x)

    public:
        Poisson(MathSolver &math);
        Eigen::VectorXd solver(std::vector<double> &S, std::vector<double> &F) override;
        double D(double x) override;
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

        double dB_du(double x, double u) override;
        double dC_du(double x, double u) override;
        double dD_du(double x, double u) override;


        double B(double x) override;
        double C(double x) override;
        double D(double x) override;
        double E(double x) override;
        double F(double x) override;

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

        double dC_du(double x, double u) override;
        double dD_du(double x, double u) override;

        double C(double x) override;
        double D(double x) override;
        double E(double x) override;
        double F(double x) override;

};

#endif
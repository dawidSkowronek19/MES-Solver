#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <string>
#include "../mesh/Point.hpp"


class Math{};

class Physics{

    public:
        Physics(Math &m_math);
        virtual ~Physics() = default;
        virtual double A(Position r);
        virtual double B(Position r);
        virtual double C(Position r);
        virtual double D(Position r);
        virtual double E(Position r);

        virtual Position v(Position r);

        virtual Eigen::VectorXd stationary_solver(const Eigen::MatrixXd& S, const Eigen::VectorXd& F);
        virtual std::pair<Eigen::VectorXd, Eigen::MatrixXd> eigen_solver(const Eigen::MatrixXd& S, const Eigen::MatrixXd& M);
    protected:
        Math &m_math;

    private:
        const double *m_t = nullptr;

};

class Poisson : public Physics{
    public:
        Poisson(Math &math);
        Eigen::VectorXd stationary_solver(const Eigen::MatrixXd& S, const Eigen::VectorXd& F) override;
        double A(Position r) override;
        double C(Position r) override;
};

class GeneralPDE : public Physics{
    public:
        GeneralPDE(Math &math);
        Eigen::VectorXd stationary_solver(const Eigen::MatrixXd& S, const Eigen::VectorXd& F) override;
        std::pair<Eigen::VectorXd, Eigen::MatrixXd> eigen_solver(const Eigen::MatrixXd& S, const Eigen::MatrixXd& M) override;

        double A(Position r) override;
        double B(Position r) override;
        double C(Position r) override;
        double D(Position r) override;
        double E(Position r) override;
        
        Position v(Position r) override;
};

class GeneralSymetricPDE : public Physics{

};
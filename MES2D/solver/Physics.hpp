#ifndef PHYSICS_HPP
#define PHYSICS_HPP

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <string>
#include <functional>
#include "../mesh/Point.hpp"
#include "../mesh/Grid.hpp"

//parents

class BilinearOperator{

    public:
        BilinearOperator(const std::vector<Triangle> &triangles, const std::vector<Position> &m_points);
        virtual ~BilinearOperator() = default;
        virtual Eigen::MatrixXd localMatrix() = 0;

    protected:
        const std::vector<Triangle> &m_triangles;
        const std::vector<Position> &m_points;
};

class LinearOperator{

    public:
        LinearOperator(const std::vector<Triangle> &triangles, const std::vector<Position> &m_points);
        virtual ~LinearOperator() = default;
        virtual Eigen::VectorXd localVector() = 0;

    protected:
        const std::vector<Triangle> &m_triangles;
        const std::vector<Position> &m_points;
};


// ==================

//childrens

class LaplaceIntegrator : public BilinearOperator{
    
    public:
        LaplaceIntegrator(const std::vector<Triangle> &Triangles, const std::vector<Position> &points,
                            const std::function<double(double)> &A);
        LaplaceIntegrator(std::function<double(double)> &A);
        Eigen::MatrixXd localMatrix() override;
    private:
        std::function<double(double)> m_A;

};

class SourceIntegrator : public LinearOperator{
    public:
        SourceIntegrator(const std::vector<Triangle> &Triangles, const std::vector<Position> &points,
                        const std::function<double(double)> &A);
        Eigen::VectorXd localVector() override;
    private:
        std::function<double(double)> m_A;
};

#endif
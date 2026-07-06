#ifndef PHYSICS_HPP
#define PHYSICS_HPP

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <string>
#include <functional>
#include "../mesh/Point.hpp"
#include "../mesh/Grid.hpp"
#include "./shapeFunction.hpp"

//parents

double Integrator();

class BilinearOperator{

    public:
        BilinearOperator(const Grid2D &Grid, const ShapeFunction &sh_func);
        virtual ~BilinearOperator() = default;
        virtual Eigen::MatrixXd S_loc(size_t element_idx) = 0;

    protected:
        const Grid2D &m_Grid;
        const ShapeFunction &m_shfunc;
};

class LinearOperator{

    public:
        LinearOperator(const Grid2D &Grid, const ShapeFunction &sh_func);
        virtual ~LinearOperator() = default;
        virtual Eigen::VectorXd F_loc(size_t element_idx) = 0;

    protected:
        const Grid2D &m_Grid;
        const ShapeFunction &m_shfunc;
};


// ==================

//childrens

class LaplaceIntegrator : public BilinearOperator{
    
    public:
        LaplaceIntegrator(const Grid2D &Grid, const ShapeFunction &sh_func  ,const std::function<std::vector<double>(Position)> &A);
        Eigen::MatrixXd S_loc(size_t element_idx) override;
    private:
        std::function<std::vector<double>(Position)> m_A;

};

class AdvectionIntegrator : public BilinearOperator{
    public: 
        AdvectionIntegrator(const Grid2D &Grid, const ShapeFunction &sh_func ,const std::function<Position(Position)> &v);
        Eigen::MatrixXd S_loc(size_t element_idx) override;
    private:
        std::function<Position(Position)> m_v;
};

class ReactionIntegrator : public BilinearOperator{
    public: 
        ReactionIntegrator(const Grid2D &Grid, const ShapeFunction &sh_func, const std::function<double(Position)> &A);
        Eigen::MatrixXd S_loc(size_t element_idx) override;
    private:
        std::function<double(Position)> m_A;
};

class SourceIntegrator : public LinearOperator{
    public:
        SourceIntegrator(const Grid2D &Grid, const ShapeFunction &sh_func, const std::function<double(Position)> &A);
        Eigen::VectorXd F_loc(size_t element_idx) override;
    private:
        std::function<double(Position)> m_A;
};



#endif
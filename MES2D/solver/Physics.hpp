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
#include "./Quadrature.hpp"

//parents

double Integrator();

class BilinearOperator{

    public:
        BilinearOperator(const ShapeFunction &sh_func);
        virtual ~BilinearOperator() = default;
        virtual void S_loc(const ElementGeometry &Geometry,  const Quadrature &quad) = 0;
        void S_clear();
        const Eigen::MatrixXd& get_S() const;

    protected:
        const ShapeFunction &m_shfunc;
        Eigen::MatrixXd m_S;
};

class LinearOperator{

    public:
        LinearOperator(const ShapeFunction &sh_func);
        virtual ~LinearOperator() = default;
        virtual void F_loc(const ElementGeometry &Geometry, const Quadrature &quad) = 0;
        void F_clear();
        const Eigen::VectorXd& get_F() const;

    protected:
        const ShapeFunction &m_shfunc;
        Eigen::VectorXd m_F;
};

class SimpleIntegrator : public BilinearOperator{
    public:
        SimpleIntegrator(const ShapeFunction &sh_func, const std::function<double(Position)> &A, const double multi);
        void S_loc(const ElementGeometry &Geometry, const Quadrature &quad) override;
    protected:
        std::function<double(Position)> m_A;
        double m_multi;
};


// ==================

//childrens

class LaplaceIntegrator : public BilinearOperator{
    
    public:
        LaplaceIntegrator(const ShapeFunction &sh_func  ,const std::function<Eigen::MatrixXd(Position)> &A);
        void S_loc(const ElementGeometry &Geometry, const Quadrature &quad) override;
    private:
        std::function<Eigen::MatrixXd(Position)> m_A;

};

class AdvectionIntegrator : public BilinearOperator{
    public: 
        AdvectionIntegrator(const ShapeFunction &sh_func ,const std::function<Eigen::Vector2d(Position)> &v);
        void S_loc(const ElementGeometry &Geometry, const Quadrature &quad) override;
    private:
        std::function<Eigen::Vector2d(Position)> m_v;
};



class SourceIntegrator : public LinearOperator{
    public:
        SourceIntegrator(const ShapeFunction &sh_func, const std::function<double(Position)> &A);
        void F_loc(const ElementGeometry &Geometry, const Quadrature &quad) override;
    private:
        std::function<double(Position)> m_A;
};


class ReactionIntegrator : public SimpleIntegrator{
    public: 
        ReactionIntegrator(const ShapeFunction &sh_func, const std::function<double(Position)> &A);
};

class MassDampingIntegrator : public SimpleIntegrator{
    public:
        MassDampingIntegrator(const ShapeFunction &sh_func, const std::function<double(Position r)> &A);
};



#endif
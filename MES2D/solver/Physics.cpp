#include "Physics.hpp"


BilinearOperator::BilinearOperator(const Grid2D &Grid, const ShapeFunction &sh_func) : m_Grid(Grid), m_shfunc(sh_func) {}
LinearOperator::LinearOperator(const Grid2D &Grid, const ShapeFunction &sh_func) : m_Grid(Grid), m_shfunc(sh_func) {}

                                
//childrens
//left

LaplaceIntegrator::LaplaceIntegrator(const Grid2D &Grid, const ShapeFunction &sh_func, const std::function<std::vector<double>(Position)> &A) : BilinearOperator(Grid, sh_func), 
                                      m_A(A) {}

AdvectionIntegrator::AdvectionIntegrator(const Grid2D &Grid, const ShapeFunction &sh_func, const std::function<Position(Position)> &v) : BilinearOperator(Grid, sh_func),
                                        m_v(v) {}

ReactionIntegrator::ReactionIntegrator(const Grid2D &Grid, const ShapeFunction &sh_func, const std::function<double(Position)> &A) : BilinearOperator(Grid, sh_func),
                                        m_A(A) {}


//right
SourceIntegrator::SourceIntegrator(const Grid2D &Grid, const ShapeFunction &sh_func, const std::function<double(Position)> &A) : LinearOperator(Grid, sh_func), 
                                    m_A(A) {}
                            



// ==================== LOCAL MATRIX =================

// Laplace

Eigen::MatrixXd LaplaceIntegrator::S_loc(const Jacobi &jacobi_package, const Quadrature &quad)
{
    
    const auto [J, J_inv, J_det] = jacobi_package.JacobiEssentials();
    int p = m_shfunc.get_p();
    int sh_nb = (p+1)*(p+2)/2;

    Eigen::MatrixXd S;

    for (int i=0; i<sh_nb; i++)
    {
        for (int j=0; j<sh_nb; j++)
        {
            for (int q=0; q<quad.get_point_numb(); q++)
            {
                double ksi = quad.get_ksi(q);
                double eta = quad.get_eta(q);

                S(i,j)+=quad.get_weight(q);
            }
        }
    }
}
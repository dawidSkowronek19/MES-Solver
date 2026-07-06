#include "Physics.hpp"


BilinearOperator::BilinearOperator(const Grid2D &Grid, const ShapeFunction &sh_func) : m_Grid(Grid), m_shfunc(sh_func) {}
LinearOperator::LinearOperator(const Grid2D &Grid, const ShapeFunction &sh_func) : m_Grid(Grid), m_shfunc(sh_func) {}

                                
//childrens
//left

LaplaceIntegrator::LaplaceIntegrator(const Grid2D &Grid, const ShapeFunction &sh_func, const std::function<double(Position)> &A) : BilinearOperator(Grid, sh_func), 
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

Eigen::MatrixXd LaplaceIntegrator::S_loc(size_t element_idx)
{

    
}
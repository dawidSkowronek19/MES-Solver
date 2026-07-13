#include "Physics.hpp"

//Parents
BilinearOperator::BilinearOperator(const ShapeFunction &sh_func) :  m_shfunc(sh_func) {
    int p = m_shfunc.get_p();
    int sh_nb = (p+1)*(p+2)/2;

    m_S = Eigen::MatrixXd::Zero(sh_nb, sh_nb);
}
LinearOperator::LinearOperator(const ShapeFunction &sh_func) : m_shfunc(sh_func) {
    int p = m_shfunc.get_p();
    int sh_nb = (p+1)*(p+2)/2;

    m_F = Eigen::VectorXd::Zero(sh_nb);
}

void BilinearOperator::S_clear(){m_S.setZero();}
void LinearOperator::F_clear(){m_F.setZero();}
const Eigen::MatrixXd& BilinearOperator::get_S() const {return m_S;}
const Eigen::VectorXd& LinearOperator::get_F() const {return m_F;}

                                
//childrens
//left

LaplaceIntegrator::LaplaceIntegrator(const ShapeFunction &sh_func, const std::function<Eigen::MatrixXd(Position)> &A) : BilinearOperator(sh_func), 
                                      m_A(A) {}

AdvectionIntegrator::AdvectionIntegrator(const ShapeFunction &sh_func, const std::function<Eigen::Vector2d(Position)> &v) : BilinearOperator(sh_func),
                                        m_v(v) {}


SimpleIntegrator::SimpleIntegrator(const ShapeFunction &sh_func, const std::function<double(Position)> &A, const double multi) : BilinearOperator(sh_func),
                                    m_A(A), m_multi(multi) {}


ReactionIntegrator::ReactionIntegrator(const ShapeFunction &sh_func, const std::function<double(Position)> &A) : SimpleIntegrator(sh_func, A, -1.0) {}
MassDampingIntegrator::MassDampingIntegrator(const ShapeFunction &sh_func, const std::function<double(Position)> &A) : SimpleIntegrator(sh_func, A, 1.0) {}




//right
SourceIntegrator::SourceIntegrator(const ShapeFunction &sh_func, const std::function<double(Position)> &A) : LinearOperator(sh_func), 
                                    m_A(A) {}
                            


// ==================== LOCAL MATRIX =================

// Laplace

void LaplaceIntegrator::S_loc(const ElementGeometry &Geometry, const Quadrature &quad) // Tensor A MUST BE SYMMETRIC
{
    
    const auto [J, J_inv, J_det] = Geometry.get_JacobiEssentials();
    int p = m_shfunc.get_p();
    int sh_nb = (p+1)*(p+2)/2;
    Eigen::MatrixXd A_eff;

    for (int q=0; q<quad.get_point_numb(); q++)
    {
        double ksi = quad.get_ksi(q);
        double eta = quad.get_eta(q);
        Position r = Geometry.cartes(ksi, eta);
        A_eff=J_inv.transpose() * m_A(r) * J_inv;

        for (int i=0; i<sh_nb; i++)
        {
            for (int j=i; j<sh_nb; j++)
            {
                double dphi_dksi_I = m_shfunc.get_dphi_dksi(q, i);
                double dphi_deta_I = m_shfunc.get_dphi_deta(q,i);

                double s=J_det*quad.get_weight(q)*(
                    m_shfunc.get_dphi_dksi(q,j)*(A_eff(0,0)*dphi_dksi_I + A_eff(1,0)*dphi_deta_I) +
                    m_shfunc.get_dphi_deta(q,j)*(A_eff(0,1)*dphi_dksi_I + A_eff(1,1)*dphi_deta_I)
                );

                m_S(i,j)+=s;
                if (i!=j) m_S(j,i)+=s;
            }
        }

    }

}

void AdvectionIntegrator::S_loc(const ElementGeometry &Geometry, const Quadrature &quad)
{
    const auto [J, J_inv, J_det] = Geometry.get_JacobiEssentials();
    int p = m_shfunc.get_p();
    int sh_nb = (p+1)*(p+2)/2;
    Eigen::Vector2d v_eff;


    for (int q=0; q<quad.get_point_numb(); q++)
    {
        double ksi = quad.get_ksi(q);
        double eta = quad.get_eta(q);

        Position r = Geometry.cartes(ksi, eta);
        v_eff = m_v(r).transpose()*J_inv.transpose();

        for (int i=0; i<sh_nb; i++)
        {
            for (int j=0; j<sh_nb; j++)
            {
                m_S(i,j) -= J_det*quad.get_weight(q)*m_shfunc.get_phi(q, i) * 
                ( v_eff(0)*m_shfunc.get_dphi_dksi(q,j) + v_eff(1)*m_shfunc.get_dphi_deta(q,j) );
            }
        }

    }
}

void SimpleIntegrator::S_loc(const ElementGeometry &Geometry, const Quadrature &quad)
{
    const auto [J, J_inv, J_det] = Geometry.get_JacobiEssentials();
    int p = m_shfunc.get_p();
    int sh_nb = (p+1)*(p+2)/2;

    for (int q=0; q<quad.get_point_numb(); q++)
    {
        double ksi = quad.get_ksi(q);
        double eta = quad.get_eta(q);

        Position r = Geometry.cartes(ksi, eta);

        for (int i=0; i<sh_nb; i++)
        {
            for (int j=i; j<sh_nb; j++)
            {
                double s = m_multi*J_det*quad.get_weight(q)*m_A(r)*m_shfunc.get_phi(q,i)*m_shfunc.get_phi(q,j);
                m_S(i,j)+=s;
                if (i!=j)   m_S(j,i)+=s; 
            }
        }
    }
}

void SourceIntegrator::F_loc(const ElementGeometry &Geometry, const Quadrature &quad)
{
    const auto [J, J_inv, J_det] = Geometry.get_JacobiEssentials();
    int p = m_shfunc.get_p();
    int sh_nb = (p+1)*(p+2)/2;

    for (int q=0; q<quad.get_point_numb(); q++)
    {
        double ksi = quad.get_ksi(q);
        double eta = quad.get_eta(q);

        Position r = Geometry.cartes(ksi, eta);

        for (int i=0; i<sh_nb; i++)
        {
            m_F(i) += quad.get_weight(q) * J_det * m_A(r) * m_shfunc.get_phi(q,i);
        }

    }
}


#include "timeDep.hpp"

Newmark_beta::Newmark_beta(const Eigen::SparseMatrix<double> &S, const Eigen::SparseMatrix<double> &M, 
            const Eigen::SparseMatrix<double> &D, const Eigen::VectorXd &F, double beta, double gamma, double dt) : 
            m_S(S), m_D(D), m_M(M), m_F(F)
{
    m_beta=beta;
    m_gamma=gamma;
    m_dt=dt;
    m_size = m_F.rows();
    m_aCst.resize(6);
    m_aCst[0] = 1.0/(m_beta*m_dt*m_dt);
    m_aCst[1] = m_gamma/(m_beta*m_dt);
    m_aCst[2] = 1.0/(m_beta*m_dt);
    m_aCst[3] = 1.0/(2.0*m_beta) - 1.0;
    m_aCst[4] = m_gamma/m_beta - 1.0;
    m_aCst[5] = 0.5*m_dt*(m_gamma/m_beta-2.0);

    m_new.C.setZero(m_size);
    m_new.dC.setZero(m_size);
    m_new.d2C.setZero(m_size);

    m_old.C.setZero(m_size);
    m_old.dC.setZero(m_size);
    m_old.d2C.setZero(m_size);
}

void Newmark_beta::set_Seffective()
{

    S_effective = m_S;
    S_effective += m_aCst[0] * m_M;
    S_effective += m_aCst[1] * m_D;

    m_solver.compute(S_effective);
}

void Newmark_beta::set_initial(const Eigen::VectorXd &C, const Eigen::VectorXd &dC)
{
    m_old.C = C;
    m_old.dC = dC;
    
    Eigen::VectorXd Res = m_F;
    Res.noalias() -=m_S*m_old.C;
    Res.noalias() -=m_D*m_old.dC;

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(m_M);
    m_old.d2C = solver.solve(Res);
}


void Newmark_beta::set_Feffective()
{
    Eigen::VectorXd termM = m_aCst[0]*m_old.C + m_aCst[2]*m_old.dC + m_aCst[3]*m_old.d2C;
    Eigen::VectorXd termD = m_aCst[1]*m_old.C + m_aCst[4]*m_old.dC + m_aCst[5]*m_old.d2C;

    F_effective = m_F;
    F_effective.noalias() += m_M * termM;
    F_effective.noalias() += m_D * termD;
}

void Newmark_beta::solve()
{
    set_Feffective();
    std::swap(m_new, m_old);
    m_new.C = m_solver.solve(F_effective);
    set_Cdiv_new();
}

void Newmark_beta::set_Cdiv_new()
{
    m_new.d2C = m_aCst[0]*(m_new.C - m_old.C) - m_aCst[2]*m_old.dC - m_aCst[3]*m_old.d2C;
    m_new.dC = m_old.dC + m_dt*((1.0-m_gamma)*m_old.d2C + m_gamma*m_new.d2C);
}

const C_vec& Newmark_beta::get_C() const {return m_new;}
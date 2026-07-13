#include "timeDep.hpp"

Newmark_beta::Newmark_beta(const Eigen::SparseMatrix<double> &S, const Eigen::SparseMatrix<double> &M, 
            const Eigen::SparseMatrix<double> &D, const Eigen::VectorXd &F) : m_S(S), m_M(M), m_D(D), m_F(F) {}

void Newmark_beta::setC(Eigen::VectorXd &C)
{
    std::swap(m_Cnew, m_Cold);
    std::swap(m_Cnew, C);
}

void Newmark_beta::setGamma(const double gamma) {m_gamma=gamma;}
void Newmark_beta::setBeta(const double beta) {m_beta=beta;}
void Newmark_beta::setDt(const double dt) {m_dt=dt;}



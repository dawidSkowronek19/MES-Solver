#ifndef TIMEDEP_HPP
#define TIMEDEP_HPP
#include <Eigen/Sparse>

struct C_vec{
    Eigen::VectorXd C, dC, d2C;
};


class Newmark_beta{
    public:
        Newmark_beta(const Eigen::SparseMatrix<double> &m_S, const Eigen::SparseMatrix<double> &m_M, 
            const Eigen::SparseMatrix<double> &m_D, const Eigen::VectorXd &m_F);
        
        void setC(Eigen::VectorXd &C);
        void setGamma(const double gamma);
        void setBeta(const double gamma);
        void setDt(const double dt);
    private:
        const Eigen::SparseMatrix<double> m_S;
        const Eigen::SparseMatrix<double> m_M;
        const Eigen::SparseMatrix<double> m_D;
        const Eigen::VectorXd m_F;
        C_vec m_new;
        C_vec m_old;

        double m_beta;
        double m_gamma;
        double m_dt;
        
        
};


#endif
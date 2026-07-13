#ifndef TIMEDEP_HPP
#define TIMEDEP_HPP
#include <Eigen/Sparse>

struct C_vec{
    Eigen::VectorXd C, dC, d2C;
};


class Newmark_beta{
    public:
        Newmark_beta(const Eigen::SparseMatrix<double> &S, const Eigen::SparseMatrix<double> &M, 
            const Eigen::SparseMatrix<double> &D, const Eigen::VectorXd &m_F, double beta, double gamma, double dt);
        void set_Seffective();
        void set_initial(const Eigen::VectorXd &C, const Eigen::VectorXd &dC);
        void solve();
        const C_vec& get_C() const;
        


        private:
        std::vector<double> m_aCst;
        const Eigen::SparseMatrix<double> &m_S;
        const Eigen::SparseMatrix<double> &m_D;
        const Eigen::SparseMatrix<double> &m_M;
        const Eigen::VectorXd &m_F;

        Eigen::VectorXd F_effective;
        Eigen::SparseMatrix<double> S_effective;
        Eigen::SparseLU<Eigen::SparseMatrix<double>> m_solver;

        C_vec m_new;
        C_vec m_old;

        double m_beta;
        double m_gamma;
        double m_dt;
        int m_size;



        
        void set_Feffective();
        void set_Cdiv_new();
        
};


#endif
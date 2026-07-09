#include "Assembler.hpp"

Assembler::Assembler(int size)
{
    m_F = Eigen::VectorXd::Zero(size);
}



void Assembler::add_Smatrix(const Eigen::MatrixXd& S_loc, const std::vector<int> &dof_indices)
{
    int localSize = S_loc.rows();
    for (int i=0; i<localSize; i++)
    {
        int global_I = dof_indices[i];
        for (int j=0; j<localSize; j++)
        {
            int global_J = dof_indices[j];
            m_triplets.push_back(Eigen::Triplet<double>(global_I, global_J, S_loc(i, j)));
        }
    }
}

void Assembler::add_Fload(const Eigen::VectorXd& F_loc, const std::vector<int> &dof_indices)
{
    int localSize = F_loc.rows();

    for (int i=0; i<localSize; i++)
    {
        int global_I = dof_indices[i];
        m_F[global_I] = F_loc[i];
    }
}

void Assembler::apply_Dirichlet_BC(const std::map<int, double>& boundary_values)
{
    const double penalty = 1e30;

    for (const auto& [dof_idx, dirichlet_val] : boundary_values)
    {
        m_triplets.push_back(Eigen::Triplet<double>(dof_idx, dof_idx, penalty));
        m_F(dof_idx) = dirichlet_val * penalty;
    }
}

void Assembler::assemble(const int total_dofs)
{
    m_S = Eigen::SparseMatrix<double>(total_dofs, total_dofs);
    m_S.setFromTriplets(m_triplets.begin(), m_triplets.end());
}


const Eigen::SparseMatrix<double>& Assembler::get_stiffnessM() const {return m_S;}
const Eigen::VectorXd& Assembler::get_loadV() const {return m_F;}
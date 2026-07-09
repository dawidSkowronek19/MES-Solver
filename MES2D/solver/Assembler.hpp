#ifndef ASSEMBLER_HPP
#define ASSEMBLER_HPP
#include <vector>
#include <map>
#include <Eigen/Sparse>
#include "./Physics.hpp"
#include "../mesh/Grid.hpp"
#include "../mesh/DoFHandler.hpp"


class Assembler{
    public:
        Assembler();
        void add_Smatrix(const Eigen::MatrixXd& S_loc, const std::vector<int> &dof_indices);
        void add_Fload(const Eigen::VectorXd& F_loc, const std::vector<int> &dof_indices);
        void assemble(const int total_dofs);
        void apply_Dirichlet_BC(const std::map<int, double>& boundary_values);
        //void apply_Neuman_BC(const std::vector<int> &boundary_nodes_idx);

        const Eigen::SparseMatrix<double>& get_stiffnessM() const;
        const Eigen::VectorXd& get_loadV() const;

    private:

        Eigen::SparseMatrix<double> m_S;
        Eigen::VectorXd m_F;
        std::vector<Eigen::Triplet<double>> m_triplets;


};

#endif
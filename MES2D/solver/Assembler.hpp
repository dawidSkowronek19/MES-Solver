#ifndef ASSEMBLER_HPP
#define ASSEMBLER_HPP
#include <vector>
#include <memory>
#include <Eigen/Sparse>
#include "./Physics.hpp"
#include "../mesh/Grid.hpp"

class Assembler{
    public:
        Assembler(const Grid2D &grid);
        void add_stiffness_term(std::shared_ptr<BilinearOperator> op);
        void add_load_term(std::shared_ptr<LinearOperator> op);
        void assemble();
        void apply_Dirichlet_BC(const std::vector<int> &boundary_nodes_idx);
        void apply_Neuman_BC(const std::vector<int> &boundary_nodes_idx);

        const Eigen::SparseMatrix<double>& get_stiffnessM() const;
        const Eigen::VectorXd& get_loadV() const;

    private:
        const Grid2D &m_grid;
        std::vector<std::shared_ptr<BilinearOperator>> m_bilinearOP;
        std::vector<std::shared_ptr<LinearOperator>> m_linearOP;

        Eigen::SparseMatrix<double> m_S;
        Eigen::VectorXd m_F;
        std::vector<Eigen::Triplet<double>> m_triplets;


};

#endif
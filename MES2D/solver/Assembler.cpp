#include "Assembler.hpp"

Assembler::Assembler(const Grid2D& grid) : m_grid(grid) {}

void Assembler::add_stiffness_term(std::shared_ptr<BilinearOperator> op)
{
    m_bilinearOP.push_back(op);
}
void Assembler::add_load_term(std::shared_ptr<LinearOperator> op)
{
    m_linearOP.push_back(op);
}

void Assembler::assemble()
{
    
}
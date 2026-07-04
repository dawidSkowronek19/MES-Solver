#include "Physics.hpp"


BilinearOperator::BilinearOperator(const std::vector<Triangle> &Triangles, const std::vector<Position> &points) :
                                    m_triangles(Triangles), m_points(points) {}
LinearOperator::LinearOperator(const std::vector<Triangle> &Triangles, const std::vector<Position> &points) :
                                    m_triangles(Triangles), m_points(points) {}

                                


LaplaceIntegrator::LaplaceIntegrator(const std::vector<Triangle> &Triangles, const std::vector<Position> &points,
                                      const std::function<double(double)> &A) : BilinearOperator(Triangles, points), 
                                      m_A(A) {}


SourceIntegrator::SourceIntegrator(const std::vector<Triangle> &Triangles, const std::vector<Position> &points,
                                    const std::function<double(double)> &A) : LinearOperator(Triangles, points), 
                                    m_A(A) {}
                            


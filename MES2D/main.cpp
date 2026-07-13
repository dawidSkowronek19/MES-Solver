#include "./mesh/Grid.hpp"
#include "./solver/shapeFunction.hpp"
#include "./solver/Quadrature.hpp"
#include "./solver/Physics.hpp"
#include "./solver/Assembler.hpp"
#include "./solver/Quadrature.hpp"

#include <memory>
#include <fstream>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

Eigen::Matrix2d laplace_tensor(Position r)
{
    Eigen::Matrix2d A;
    A<<1.0, 0.0,
        0.0, 1.0;
    return A;
}
double load_scalar(Position r)
{
    if (r.x()*r.x() + r.y()*r.y() < 0.04) return 50.0;
    return 0.0;
}

double bc_func(Position r)
{
    return 0.0;
}

Eigen::Vector2d v(Position r)
{
    return {r.y(), -r.x()};
}

int main()
{
    GeometryParameters geoparam{
        200,
        2*sin(M_PI/200),
        2*0.5*sin(M_PI/200),
        100
    };

    int p = 1;

    Eigen::MatrixXd S_local;
    Eigen::MatrixXd F_local;

    Grid2D circle(geoparam);
    circle.triangular();
    circle.create_neighboursList();
    circle.relaxGrid();

    circle.saveTrianglesPoints();
    circle.savePointsList();

    std::cout<<"quadrature\n";
    Quadrature quad = QuadratureFactory::refTriangle(p);
    std::cout<<"shape functions\n";
    ShapeFunction sh_func(p);
    sh_func.set_cached_values(quad.get_integrationPoints());
    std::cout<<"degree of freedom\n";
    DoFHandler dof(circle, p);
    dof.countNodes();
    dof.set_boundary_dofs(geoparam.EdgeNodesNumber, bc_func);
    std::cout << "DOF: " << dof.get_totalDOF() << std::endl;
    auto bc_map = dof.get_boundary_dofs();

    std::cout<<"operators\n";
    std::vector<std::shared_ptr<BilinearOperator>> bilinear_ops;
    std::vector<std::shared_ptr<LinearOperator>> linear_ops;

    bilinear_ops.push_back(std::make_shared<LaplaceIntegrator>(sh_func, laplace_tensor));
    //bilinear_ops.push_back(std::make_shared<AdvectionIntegrator>(sh_func, v));


    linear_ops.push_back(std::make_shared<SourceIntegrator>(sh_func, load_scalar));
    std::cout<<"assembler init\n";
    Assembler assembler(dof.get_totalDOF());
    const std::vector<std::vector<int>>& elements_dofs = dof.get_nodesID();
    std::cout<<"local element loop\n";
    for (int el_idx=0; el_idx<circle.getElementNb(); el_idx++)
    {
        Triangle triangle = circle.getTriangle(el_idx);
        ElementPointPositions acc_el_points;
        acc_el_points.p1 = circle.getPoint(triangle.p1);
        acc_el_points.p2 = circle.getPoint(triangle.p2);
        acc_el_points.p3 = circle.getPoint(triangle.p3);

        ElementGeometry geometry(acc_el_points);
        //std::cout<<"Jacobi\n";
        geometry.calcJacobi();
        //std::cout<<"local operators\n";
        for (auto& op : bilinear_ops)
        {
            op->S_clear();
            op->S_loc(geometry, quad);

            assembler.add_Smatrix(op->get_S(), elements_dofs[el_idx]);
        }

        for (auto& op : linear_ops)
        {
            op->F_clear();
            op->F_loc(geometry, quad);

            assembler.add_Fload(op->get_F(), elements_dofs[el_idx]);
        }
    }
    std::cout<<"boundary conditions\n";
    assembler.apply_Dirichlet_BC(bc_map);
    std::cout<<"assembling\n";
    assembler.assemble(dof.get_totalDOF());
    std::cout<<"solving\n";

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(assembler.get_stiffnessM());

    if(solver.info() != Eigen::Success) {
        std::cerr << "Decomposition failed!" << std::endl;
        return 1;
    }

    Eigen::VectorXd C = solver.solve(assembler.get_loadV());

    if(solver.info() != Eigen::Success) {
        std::cerr << "Solving failed!" << std::endl;
        return 1;
    }

    //for test
    std::ofstream file_C("./outdir/C.dat");
    int vis_resolution = 10; 

    for (int el_idx = 0; el_idx < circle.getElementNb(); el_idx++)
    {
        Triangle triangle = circle.getTriangle(el_idx);
        ElementPointPositions acc_el_points;
        acc_el_points.p1 = circle.getPoint(triangle.p1);
        acc_el_points.p2 = circle.getPoint(triangle.p2);
        acc_el_points.p3 = circle.getPoint(triangle.p3);

        ElementGeometry geometry(acc_el_points);
        const std::vector<int>& local_dofs = elements_dofs[el_idx];

        for (int i = 0; i <= vis_resolution; i++)
        {
            for (int j = 0; j <= vis_resolution - i; j++)
            {
                double ksi = static_cast<double>(i) / vis_resolution;
                double eta = static_cast<double>(j) / vis_resolution;

                Position r = geometry.cartes(ksi, eta);
                double u_val = 0.0;
                for (int d = 0; d < local_dofs.size(); d++)
                {
                    int global_dof = local_dofs[d];
                    u_val += C(global_dof) * sh_func.phi(ksi, eta, d);
                }
                file_C << r.x() << " " << r.y() << " " << u_val << "\n";
            }
        }
    }   

    file_C.close();
    std::cout << "done!\n"; 

}   

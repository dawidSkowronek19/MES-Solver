#include "./mesh/Grid.hpp"
#include "./solver/shapeFunction.hpp"
#include "./solver/Quadrature.hpp"
#include "./solver/Physics.hpp"
#include "./solver/Assembler.hpp"
#include "./solver/Quadrature.hpp"
#include "./solver/timeDep.hpp"

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
    //if (r.x()*r.x() + r.y()*r.y() < 0.04) return 50.0;
    return 0.0;
}

double bc_func(Position r)
{
    return 1.0;
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
    std::vector<std::shared_ptr<BilinearOperator>> mass_ops;
    std::vector<std::shared_ptr<BilinearOperator>> damp_ops;
    std::vector<std::shared_ptr<LinearOperator>> linear_ops;

    bilinear_ops.push_back(std::make_shared<LaplaceIntegrator>(sh_func, laplace_tensor));
    linear_ops.push_back(std::make_shared<SourceIntegrator>(sh_func, load_scalar));
    mass_ops.push_back(std::make_shared<MassDampingIntegrator>(sh_func, [](Position r){return 1.0;} ));
    damp_ops.push_back(std::make_shared<MassDampingIntegrator>(sh_func, [](Position r){return 0.0;} ));

    std::cout<<"Assembler init\n";
    Assembler main_assembler(dof.get_totalDOF());
    Assembler mass_assembler(dof.get_totalDOF());
    Assembler damp_assembler(dof.get_totalDOF());


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

            main_assembler.add_Smatrix(op->get_S(), elements_dofs[el_idx]);
        }

        for (auto& op : linear_ops)
        {
            op->F_clear();
            op->F_loc(geometry, quad);

            main_assembler.add_Fload(op->get_F(), elements_dofs[el_idx]);
        }
        
        for (auto& op : mass_ops)
        {
            op->S_clear();
            op->S_loc(geometry, quad);

            mass_assembler.add_Smatrix(op->get_S(), elements_dofs[el_idx]);
        }

        for (auto& op : damp_ops)
        {
            op -> S_clear();
            op -> S_loc(geometry, quad);

            damp_assembler.add_Smatrix(op->get_S(), elements_dofs[el_idx]);
        }
    }

    std::cout<<"boundary conditions\n";
    main_assembler.apply_Dirichlet_BC(bc_map);
    mass_assembler.apply_Dirichlet_BC(bc_map);
    damp_assembler.apply_Dirichlet_BC(bc_map);
    std::cout<<"assembling\n";
    main_assembler.assemble(dof.get_totalDOF());
    mass_assembler.assemble(dof.get_totalDOF());
    damp_assembler.assemble(dof.get_totalDOF());

    double t=0.0;
    double dt = 0.01; 
    double tmax=1.0;
    int freq = 5;
    int iter=0;
    int vis_resolution = 5;
    Newmark_beta newmark(main_assembler.get_stiffnessM(), mass_assembler.get_stiffnessM(), damp_assembler.get_stiffnessM(), main_assembler.get_loadV(),
    0.25, 0.5, dt);
    newmark.set_Seffective();
    newmark.set_initial(Eigen::VectorXd::Zero(dof.get_totalDOF()), Eigen::VectorXd::Zero(dof.get_totalDOF()));
    while(t<tmax)
    {
        iter++;
        t+=dt;
        newmark.solve();
        auto [C, dC, d2C] = newmark.get_C();
        if (iter%freq==0)
        {
            std::string filename = "./outdir/u/sol_"+std::to_string(iter) + ".dat";
            std::ofstream file_C(filename);

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
        }

    }

    return 0;
}
#include "Solver.h"

//===================== SOLVER PART ==============================

//===================== CONSTRUCTOR & DESTRUCTOR =================

Solver::Solver(Grid_1D &grid, ShapeFunction &shapefunction, Physics &physics) : m_grid(grid), m_shapefunction(shapefunction), m_physics(physics)
{
    std::cout<<std::string(60, '=')<<"\n";
    std::cout<<"\n\n\n \t\t\tSOLVER SECTION\n\n\n";
    m_dxSave=0.01;
    m_integrationOrder=ceil((m_shapefunction.get_deg()*m_shapefunction.get_deg()+1)/2.0);
    

    const gsl_integration_fixed_type * T = gsl_integration_fixed_legendre;

    m_work = gsl_integration_fixed_alloc(T, m_integrationOrder, -1.0, 1.0, 0.0, 0.0);
    m_local_nodes=gsl_integration_fixed_nodes(m_work);
    m_local_weights=gsl_integration_fixed_weights(m_work);

    std::cout<<"\t\t"<<std::string(4, '#')<<" SOLVER PARAMETERS "<<std::string(4, '#')<<"\n\n";
    std::cout<<"\tdim[phi(x)] = "<<m_shapefunction.get_deg()<<"\n";
    std::cout<<"\tsaving step size = "<<m_dxSave<<"\n";
    std::cout<<"\tIntegration Order = "<<m_integrationOrder<<"\n\n";
    std::cout<<"\t\t"<<"# SOLVER PARAMETERS END #\n\n";

}

Solver::~Solver() 
{
    gsl_integration_fixed_free(m_work);
    std::cout<<"\n\n\n \t\t\tSOLVER SECTION DONE\n";
    std::cout<<std::string(60, '=')<<"\n";
}

double Solver::retOne(double x) {return 1.0;}
//===========================================================

//====================== SOLVER ============================

// 1D REAL


// ****************** LOCAL ELEMENT *******************
void Solver::local_1D_linear_StifnessMatrix(std::vector<double> &S, int m)
{
    double Jm=(m_grid.get_NodPosition(m+1)-m_grid.get_NodPosition(m))/2.0;
    int p = m_shapefunction.get_deg()+1;
    for (int j=0; j<p; j++)
    {
        for (int i=0; i<p; i++)
        {
            int l = i+p*j;

            for (int k=0; k<m_integrationOrder; k++)
            {
                double ksi = m_local_nodes[k]; 
                double x_ksi = 0.5*((m_grid.get_NodPosition(m+1)-m_grid.get_NodPosition(m))*ksi + (m_grid.get_NodPosition(m+1)+m_grid.get_NodPosition(m)));
                S[l]+=m_local_weights[k]*(-m_shapefunction.Rphi_1D_deriv(ksi, i)*m_shapefunction.Rphi_1D_deriv(ksi,j)/Jm + 
                m_physics.B(x_ksi)*m_shapefunction.Rphi_1D_deriv(ksi,i)*m_shapefunction.Rphi_1D(ksi,j)+
                m_physics.C(x_ksi)*m_shapefunction.Rphi_1D(ksi,i)*m_shapefunction.Rphi_1D(ksi,j)*Jm);
            }

        }
    }
}

void Solver::local_1D_nonlinear_StifnessMatrix(std::vector<double> &S, int m)
{
    double Jm=(m_grid.get_NodPosition(m+1)-m_grid.get_NodPosition(m))/2.0;
    int p = m_shapefunction.get_deg()+1;
    for (int j=0; j<p; j++)
    {
        for (int i=0; i<p; i++)
        {
            int l = i+p*j;

            for (int k=0; k<m_integrationOrder; k++)
            {
                double ksi = m_local_nodes[k]; 
                double x_ksi = 0.5*((m_grid.get_NodPosition(m+1)-m_grid.get_NodPosition(m))*ksi + (m_grid.get_NodPosition(m+1)+m_grid.get_NodPosition(m)));
                double u_ksi = 0.0;
                for (int k = 0; k < p; k++) 
                {
                    u_ksi += m_Q[m*(p-1)+k]*m_shapefunction.Rphi_1D(ksi, k);
                }
                S[l]+=m_local_weights[k]*(-m_shapefunction.Rphi_1D_deriv(ksi, i)*m_shapefunction.Rphi_1D_deriv(ksi,j)/Jm + 
                m_physics.B(x_ksi, u_ksi)*m_shapefunction.Rphi_1D_deriv(ksi,i)*m_shapefunction.Rphi_1D(ksi,j)+
                m_physics.C(x_ksi, u_ksi)*m_shapefunction.Rphi_1D(ksi,i)*m_shapefunction.Rphi_1D(ksi,j)*Jm);
            }

        }
    }
}

void Solver::local_1D_linear_LoadVector(std::vector<double> &F, int m)
{
    double Jm=(m_grid.get_NodPosition(m+1)-m_grid.get_NodPosition(m))/2.0;
    int p = m_shapefunction.get_deg()+1;
    for (int j=0; j<p; j++)
    {
        for (int k=0; k<m_integrationOrder; k++)
        {
            double ksi = m_local_nodes[k]; 
            double x_ksi = 0.5*((m_grid.get_NodPosition(m+1)-m_grid.get_NodPosition(m))*ksi + (m_grid.get_NodPosition(m+1)+m_grid.get_NodPosition(m)));

            F[j]+=m_local_weights[k]*Jm*m_physics.D(x_ksi)*m_shapefunction.Rphi_1D(ksi,j);
        }
    }
}



void Solver::local_1D_nonlinear_LoadVector(std::vector<double> &F, int m)
{
    double Jm=(m_grid.get_NodPosition(m+1)-m_grid.get_NodPosition(m))/2.0;
    int p = m_shapefunction.get_deg()+1;
    for (int j=0; j<p; j++)
    {
        for (int k=0; k<m_integrationOrder; k++)
        {
            double ksi = m_local_nodes[k]; 
            double x_ksi = 0.5*((m_grid.get_NodPosition(m+1)-m_grid.get_NodPosition(m))*ksi + (m_grid.get_NodPosition(m+1)+m_grid.get_NodPosition(m)));
            double u_ksi = 0.0;
            for (int k = 0; k < p; k++) {
                u_ksi += m_Q[m*(p-1)+k]*m_shapefunction.Rphi_1D(ksi, k);
            }
            F[j]+=m_local_weights[k]*Jm*m_physics.D(x_ksi, u_ksi)*m_shapefunction.Rphi_1D(ksi,j);
        }
    }
}

void Solver::local_1D_StiffnessTangentMatrix(std::vector<double> &S_T, int m)
{
    double Jm=(m_grid.get_NodPosition(m+1)-m_grid.get_NodPosition(m))/2.0;
    int p = m_shapefunction.get_deg()+1;
    for (int j=0; j<p; j++)
    {
        for (int i=0; i<p; i++)
        {
            int l = i+p*j;
            double S_basic=0.0;
            double S_extra=0.0;
            for (int int_idx=0; int_idx<m_integrationOrder; int_idx++) //integrate
            {
                double ksi=m_local_nodes[int_idx];
                double x_ksi = 0.5*((m_grid.get_NodPosition(m+1)-m_grid.get_NodPosition(m))*ksi + (m_grid.get_NodPosition(m+1)+m_grid.get_NodPosition(m)));
                double u_ksi=0.0;
                for (int k = 0; k < p; k++) 
                {
                    u_ksi += m_Q[m*(p-1)+k]*m_shapefunction.Rphi_1D(ksi, k);
                }
                S_basic+=m_local_weights[int_idx]*(-m_shapefunction.Rphi_1D_deriv(ksi, i)*m_shapefunction.Rphi_1D_deriv(ksi,j)/Jm + 
                        m_physics.B(x_ksi, u_ksi)*m_shapefunction.Rphi_1D_deriv(ksi,i)*m_shapefunction.Rphi_1D(ksi,j)+
                        (-m_physics.dD_du(x_ksi,u_ksi)+m_physics.C(x_ksi,u_ksi))*m_shapefunction.Rphi_1D(ksi,i)*m_shapefunction.Rphi_1D(ksi,j)*Jm);
            
                double sum_tmp=0.0;
                for (int k=0; k<p; k++)
                {
                    sum_tmp+=m_Q[m*m_shapefunction.get_deg()+k]*(m_physics.dB_du(x_ksi, u_ksi)*m_shapefunction.Rphi_1D_deriv(ksi,i)*m_shapefunction.Rphi_1D(ksi,j)*m_shapefunction.Rphi_1D(ksi,k)
                            + m_physics.dC_du(x_ksi,u_ksi)*m_shapefunction.Rphi_1D(ksi,i)*m_shapefunction.Rphi_1D(ksi,j)*m_shapefunction.Rphi_1D(ksi,k)*Jm);
                }
                S_extra+=m_local_weights[int_idx]*sum_tmp;
                 
            }

            S_T[l] = S_basic+S_extra;
        }

    }
    
}


void Solver::local_1D_MassMatrix(std::vector<double> &M, int m, std::function<double(double)> f)
{
    double Jm=(m_grid.get_NodPosition(m+1)-m_grid.get_NodPosition(m))/2.0;
    int p = m_shapefunction.get_deg()+1;
    for (int j=0; j<p; j++)
    {
        for (int i=0; i<p; i++)
        {
            int l = i+p*j;

            for (int k=0; k<m_integrationOrder; k++)
            {
                double ksi = m_local_nodes[k]; 
                double x_ksi = 0.5*((m_grid.get_NodPosition(m+1)-m_grid.get_NodPosition(m))*ksi + (m_grid.get_NodPosition(m+1)+m_grid.get_NodPosition(m)));
                M[l]+=m_local_weights[k]*(m_shapefunction.Rphi_1D(ksi,i)*f(x_ksi)*m_shapefunction.Rphi_1D(ksi,j)*Jm);
            }

        }
    }
}

void Solver::local_1D_MassMatrix(std::vector<double> &M, int m)
{
    local_1D_MassMatrix(M, m, retOne);
}
// ****************************************************

// ****************** BOUNDARY CONDITIONS *******************

void Solver::boundaryConditions(std::string work_type)
{
    int p = m_shapefunction.get_deg()+1;
    int N = (p-1)*m_grid.get_elementNumber()+1;

    std::cout<<"\tIncluding boundary conditions...";
    if (work_type=="stationary")
    {
        for (const auto &bc : m_grid.get_BC())
        {
            int m = bc.m;
            int I_glob = m*(p-1);
            if (bc.bct==1)
            {
                for (int j=0; j<N; j++)
                {
                    int L_glob = I_glob + j*N;
                    int L_glob_tr = j+I_glob*N;
                    m_F1D[j]-=bc.bc_value*m_S1D[L_glob];
                    m_S1D[L_glob]=0.0;
                    m_S1D[L_glob_tr]=0.0;

                    if(I_glob==j)
                        m_S1D[L_glob]=1.0;
                }

                m_F1D[I_glob]=bc.bc_value;
            }

            else if (bc.bct==2)
            {

                if (I_glob == 0) 
                {
                m_F1D[I_glob] += bc.bc_value; 
                }   
                else 
                {
                    m_F1D[I_glob] -= bc.bc_value; 
                }
            }
        }
    }

    if (work_type=="eigen")
    {
        for (const auto &bc : m_grid.get_BC())
        {
            int m = bc.m;
            int I_glob = m*(p-1);
            if (bc.bct==1)
            {
                
                for (int j=0; j<N; j++)
                {
                    int L_glob = I_glob + j*N;
                    int L_glob_tr=j+I_glob*N;

                    m_S1D[L_glob]=0.0;
                    m_M1D[L_glob]=0.0;

                    m_S1D[L_glob_tr]=0.0;
                    m_M1D[L_glob_tr]=0.0;

                    if(I_glob==j)
                    {
                        m_S1D[L_glob]=1.0;
                        m_M1D[L_glob]=1.0;
                    }
                }    


            }

            else if (bc.bct==2)
            {
                std::cout<<"NO EIGEN PROBLEM FOR NEUMAN CONDITION\n";
            }
        }
    }

    if (work_type=="nonlinear")
    {
        for (const auto &bc : m_grid.get_BC())
        {
            if (bc.bct == 1)
            {
                int p = m_shapefunction.get_deg() + 1;
                int I_glob = bc.m * (p - 1); 


                for (int j = 0; j < N; j++) 
                {
                    m_St1D[I_glob + j * N] = 0.0; 
                }

                m_St1D[I_glob + I_glob * N] = 1.0; 
                m_R1D[I_glob] = m_Q[I_glob]-bc.bc_value; 
            }
        }   
    }
    std::cout<<"DONE\n";

}


// ****************************************************

// ****************** ASSEMBLERS *******************

void Solver::Matrix_assembler(std::string work_type, std::vector<double> &matrix, std::function<double(double)> f)
{
    int p = m_shapefunction.get_deg()+1;
    std::vector<double> S_loc(p*p, 0.0);
    int N = (p-1)*m_grid.get_elementNumber()+1; 

    matrix.resize(N*N, 0.0);

    std::cout<<"\t"<<work_type<<" assemble process start... ";
    std::function<void(std::vector<double>&, int)> local_matrix_func = nullptr;

    if (work_type == "stiffness_linear") 
    {
        local_matrix_func = [this](std::vector<double>& S, int m) 
        { 
            this->local_1D_linear_StifnessMatrix(S, m); 
        };
    }
    else if (work_type == "stiffness_nonlinear") 
    {
        local_matrix_func = [this](std::vector<double>& S, int m) 
        { 
            this->local_1D_nonlinear_StifnessMatrix(S, m); 
        };
    } 
    else if (work_type == "eigen") 
    {
        local_matrix_func = [this](std::vector<double>& S, int m) 
        { 
            this->local_1D_MassMatrix(S, m); 
        };
    } 
    else if (work_type == "tangent") 
    {
        local_matrix_func = [this](std::vector<double>& S, int m) 
        { 
            this->local_1D_StiffnessTangentMatrix(S, m); 
        };
    }
    else if (work_type =="time_dependent")
    {
        local_matrix_func = [this, f](std::vector<double>& S, int m) 
        { 
            this->local_1D_MassMatrix(S, m, f); 
        };
    }
    else 
    {
        std::cout << "[ Error ]: WRONG WORK_TYPE: " << work_type << "\n";
        return; 
    }

    for (int m=0; m<m_grid.get_elementNumber(); m++)
    {
        S_loc.assign(S_loc.size(), 0.0);

    local_matrix_func(S_loc, m);

        for (int j=0; j<p; j++)
        {
            int J_glob = m*(p-1)+j;
            for(int i=0; i<p; i++)
            {
                int I_glob=m*(p-1)+i;

                int L_glob=I_glob+N*J_glob;
                int l = i+p*j;
                matrix[L_glob]+=S_loc[l];
            }
        }
    }

    std::cout<<"DONE\n";

}

void Solver::Vector_assembler(std::string work_type)
{
    int p = m_shapefunction.get_deg()+1;
    std::vector<double> F_loc(p, 0.0);
    int N = (p-1)*m_grid.get_elementNumber()+1; 
    m_F1D.resize(N, 0.0);

    std::cout<<"\tLoad assemble process start... ";
    void (Solver::*local_vector_func)(std::vector<double>&, int) = nullptr;

    if (work_type == "load_linear") 
    {
        local_vector_func = &Solver::local_1D_linear_LoadVector;
    }
    else if (work_type == "load_nonlinear") 
    {
        local_vector_func = &Solver::local_1D_nonlinear_LoadVector;
    }
        else 
    {
        std::cout << "[ Error ]: WRONG WORK_TYPE: " << work_type << "\n";
        return; 
    }

    for (int m=0; m<m_grid.get_elementNumber(); m++)
    {
        F_loc.assign(F_loc.size(), 0.0);
        (this->*local_vector_func)(F_loc, m);

        for (int j=0; j<p; j++)
        {
            int J_glob = m*(p-1)+j;
            m_F1D[J_glob]+=F_loc[j];
        }
    }

    std::cout<<"DONE\n";

}



// ****************************************************

// ****************** SOLVER STATIONARY 1 DIM *******************

void Solver::stationary_1D_linear()
{
    std::cout<<"\t\t"<<std::string(4, '#')<<" SOLVER START "<<std::string(4, '#')<<"\n\n";
    std::cout<<"\t SOLVER MODE: stationary 1D linear\n";
    
    int N = m_shapefunction.get_deg()*m_grid.get_elementNumber()+1; 
    m_Q.resize(N);
    Matrix_assembler("stiffness_linear", m_S1D);
    Vector_assembler("load_linear");
    boundaryConditions("stationary");

    Eigen::Map<Eigen::VectorXd> Q_e(m_Q.data(), N);

    Q_e = m_physics.solver(m_S1D, m_F1D);

    std::cout<<"\t\t\tPROCESS DONE\n";
}

void Solver::stationary_1D_nonlinear()
{
    std::cout<<"\t\t"<<std::string(4, '#')<<" SOLVER START "<<std::string(4, '#')<<"\n\n";
    std::cout<<"\t SOLVER MODE: stationary 1D nonlinear (Newton- Raphson method)\n";
    
    int N = m_shapefunction.get_deg()*m_grid.get_elementNumber()+1;
    int p = m_shapefunction.get_deg()+1;
    m_Q.resize(N, 1.0);
    m_St1D.resize(N*N,0.0);
    m_S1D.resize(N*N,0.0);
    m_F1D.resize(N, 0.0);
    m_R1D.resize(N,1.0);


    double epsilon = 1e-4;
    int iterrMax = 100;
    int iterr = 0;
    double error = 1.0;
    double omega = 0.2;
    Eigen::Map<Eigen::VectorXd> Q_e(m_Q.data(), N);
    Eigen::Map<Eigen::MatrixXd> St_e(m_St1D.data(), N, N);
    Eigen::Map<Eigen::MatrixXd> S_e(m_S1D.data(), N, N);
    Eigen::Map<Eigen::VectorXd> F_e(m_F1D.data(), N);
    Eigen::Map<Eigen::VectorXd> R_e(m_R1D.data(), N);

    while(true)
    {
        iterr++;
        double norm = R_e.norm();
        std::cout<<"\tnorm= "<<norm<<"\n";

        if (norm<epsilon || iterr>iterrMax)
            break;
        std::cout<<"\tIteration number = "<<iterr<<"\n";
        m_S1D.assign(m_S1D.size(), 0.0);
        m_F1D.assign(m_F1D.size(), 0.0);
        m_St1D.assign(m_St1D.size(), 0.0);
        Matrix_assembler("stiffness_nonlinear", m_S1D);
        Vector_assembler("load_nonlinear");

        R_e =S_e*Q_e-F_e;
        
        
        Matrix_assembler("tangent", m_St1D);
        boundaryConditions("nonlinear");

        Eigen::VectorXd dQ_e = St_e.partialPivLu().solve(-R_e);
        Q_e+=omega*dQ_e;

    }

    std::cout<<"DONE\n";
}

// ****************************************************

// **************** SOLVER TIME DEPENDENT 1 DIM ********************

void Solver::timeDependent_1D_linear()
{
    // WORK IN PROGRESS
    
    int N = m_shapefunction.get_deg()*m_grid.get_elementNumber()+1;
    int p = m_shapefunction.get_deg()+1;


    double beta=0.25;
    double gamma=0.5;
    double dt=0.01;

    m_Q.resize(N, 1.0);
    Eigen::Map<Eigen::VectorXd> Q_e(m_Q.data(), N);
    Eigen::Map<Eigen::VectorXd> F_e(m_F1D.data(), N);
    Eigen::VectorXd Q_e_OLD(N), d_Q_e_OLD(N), d2_Q_e_OLD(N);// zainicjalizowac warunkami poczatkowymi i brzegowymi
    std::vector<double> ME(N*N,0.0), MF(N*N,0.0);


    Eigen::Map<Eigen::MatrixXd> ME_e(ME.data(), N, N);
    Eigen::Map<Eigen::MatrixXd> MF_e(MF.data(), N, N);
    Eigen::Map<Eigen::MatrixXd> S_e(m_S1D.data(), N, N);

    double t=0.0;
    double t_max=1.0;
    int f=10;
    int time_idx=0;
    m_physics.connectTime(&t);

    while(t<=t_max)
    {
        Matrix_assembler("time_dependent", ME, [this](double x){return m_physics.E(x);});
        Matrix_assembler("time_dependent", MF, [this](double x){return m_physics.F(x);});
        Matrix_assembler("stiffness_linear", m_S1D);
        Eigen::MatrixXd S_eff_e =  1.0/(beta*dt*dt)*ME_e + (gamma/(beta*dt))*MF_e - S_e;
        Eigen::VectorXd F_eff = ((-1.0 + 1.0/(2.0*beta))*ME_e + (gamma/(2.0*beta)-1.0)*MF_e)*d2_Q_e_OLD + 
                                (1.0/(beta*dt)*ME_e-(1.0+gamma/beta)*MF_e)*d_Q_e_OLD + 
                                (1.0/(beta*dt*dt)*ME_e + gamma/(beta*dt)*MF_e)*Q_e_OLD - F_e;

        Q_e = S_eff_e.partialPivLu().solve(F_eff);

        if(time_idx%f==0)
        {
            saveSolution("outdir", "time_dependent");
        }

        //Q', Q'' computing from Newmark formulas
        for (int j=0; j<N; j++)
        {
            double d2qTMP = d2_Q_e_OLD(j);
            d2_Q_e_OLD(j) = 1.0/(beta*dt*dt)*(Q_e(j)-Q_e_OLD(j))-1.0/(beta*dt)*d_Q_e_OLD(j) + (1.0-1.0/(2.0*beta))*d2_Q_e_OLD(j);
            d_Q_e_OLD(j) = (1.0/gamma/beta)*d_Q_e_OLD(j) + (1-gamma/(2.0*beta))*dt*d2qTMP + gamma/(beta*dt)*(Q_e(j)-Q_e_OLD(j));

        }
    
    }

}

void Solver::Eigen_1D()
{
    std::cout<<"\t\t"<<std::string(4, '#')<<" SOLVER START "<<std::string(4, '#')<<"\n\n";
    std::cout<<"\t SOLVER MODE: EIGEN PROBLEM 1D\n";
    int N = m_shapefunction.get_deg()*m_grid.get_elementNumber()+1; 

    Matrix_assembler("Stiffness", m_S1D);
    Matrix_assembler("eigen", m_M1D);
    boundaryConditions("eigen");

    //std::cout<<m_S1D.size()<<" "<<m_M1D.size()<<"\n";
    auto result = m_physics.EigenSolver(m_S1D, m_M1D);

    m_eigenValues1D.resize(N);
    m_eigenVectors1D.resize(N*N);

    Eigen::Map<Eigen::VectorXd>(m_eigenValues1D.data(), N) = result.first;
    Eigen::Map<Eigen::MatrixXd>(m_eigenVectors1D.data(), N, N) = result.second;


    std::cout<<"\t\t\tPROCESS DONE\n";
}

//==============================================

//================== SOLUTION ==================

double Solver::U_1D(double x)
{
    int N = m_grid.get_elementNumber();
    int p_deg = m_shapefunction.get_deg();
    double tol=1e-12;

    for (int m = 0; m < N; m++)
    {
        double x_left = m_grid.get_NodPosition(m);
        double x_right = m_grid.get_NodPosition(m+1);

        if (x >= x_left-tol && x <= x_right+tol)
        {

            double ksi = (2.0 * (x - x_left) / (x_right - x_left)) - 1.0;
            
            double u_val = 0.0;
            for (int i = 0; i <= p_deg; i++)
            {
                int glob_idx = m * p_deg + i; 
                u_val += m_Q[glob_idx] * m_shapefunction.Rphi_1D(ksi, i);
            }
            
            return u_val;
        }
    }
    return 0.0;
}

double Solver::U_1D(double x, int mode)
{
    int N = m_grid.get_elementNumber();
    int p_deg = m_shapefunction.get_deg();
    double tol=1e-12;

    for (int m = 0; m < N; m++)
    {
        double x_left = m_grid.get_NodPosition(m);
        double x_right = m_grid.get_NodPosition(m+1);

        if (x >= x_left-tol && x <= x_right+tol)
        {
            double ksi = (2.0 * (x - x_left) / (x_right - x_left)) - 1.0;
            
            double u_val = 0.0;
            for (int i = 0; i <= p_deg; i++)
            {
                int glob_idx = m * p_deg + i;

                int l = glob_idx + mode*(p_deg*m_grid.get_elementNumber()+1); 
                u_val += m_eigenVectors1D[l] * m_shapefunction.Rphi_1D(ksi, i);
            }
            
            return u_val;
        }
    }
    return 0.0;
}

//==============================================

//================== SAVE DATA =================

void Solver::saveSolution(std::string outdir, std::string work_type)
{
    if (work_type=="stationary")
    {
        std::cout<<"\t\t\tSaving data... ";
        auto [x_start, x_end] = m_grid.get_RodInfo();
        double x=x_start;

        std::ofstream file(outdir+"/u"+std::to_string(m_physics.getCurrentTime())+".dat");

        while (x<=x_end)
        {
            file<<x<<" "<<U_1D(x)<<"\n";
            x+=m_dxSave;
        }
    
        file<<x_end<<" "<<U_1D(x_end)<<"\n";

        file.close();
        std::cout<<"DONE\n";
    }

    else if (work_type=="eigen")
    {
        
        std::ofstream file_eigenvalue(outdir+"/eigenvalue.dat");
        
        for (int mode=0; mode<m_shapefunction.get_deg()*m_grid.get_elementNumber()+1; mode++)
        {
            auto [x_start, x_end] = m_grid.get_RodInfo();
            double x=x_start;

            std::ofstream file_eigenvector(outdir+"/psi_" + std::to_string(mode) +".dat");

            while(x<x_end)
            {
                file_eigenvector<<x<<" "<<U_1D(x, mode)<<"\n";
                x+=m_dxSave;
            }
            
            file_eigenvector<<x_end<<" "<<U_1D(x_end, mode)<<"\n";
            file_eigenvector.close();

            file_eigenvalue<<mode<<" "<<m_eigenValues1D[mode]<<"\n";
        }

        file_eigenvalue.close();
    
    }

}
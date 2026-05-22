
#include "Solver.h"
//===================== SOLVER PART ==============================

//===================== CONSTRUCTOR & DESTRUCTOR =================

Solver::Solver(Grid_1D &grid) : m_grid(grid)
{
    std::cout<<"Solver start...\n";
    m_dxSave=0.01;
    m_shapeFunctionDeg=1;
    m_integrationOrder=ceil((m_shapeFunctionDeg*m_shapeFunctionDeg+1)/2.0);

    const gsl_integration_fixed_type * T = gsl_integration_fixed_legendre;
    m_work = gsl_integration_fixed_alloc(T, m_integrationOrder, -1.0, 1.0, 0.0, 0.0);
    
    m_local_nodes=gsl_integration_fixed_nodes(m_work);
    m_local_weights=gsl_integration_fixed_weights(m_work);


}

Solver::~Solver() 
{
    gsl_integration_fixed_free(m_work);
    std::cout<<"Solver end\n";
}

//================================================================

//====================== EQUATION FUNCTIONS =====================

double Solver::B(double x) //   [d2/dx +B(x) d/dx  +C(x)] u(x) =D(x)     [d2/dx +B(x) d/dx  +C(x)] u(x) =lambda u(x)
{
    return 0.0;
}
double Solver::C(double x)
{
    return 0.0;
}
double Solver::D(double x)
{
    return 1.0;
}


//===========================================================

//====================== SOLVER ============================

// 1D REAL

std::vector<double> Solver::real_thomas(std::vector<double> &a, std::vector<double> &b, std::vector<double>& c, std::vector<double> &d)
{
    /*
        Ax=d => a_i * x_{i-1} + b_i * x_i + c_i * x_{i+1} = d_i  
    */
    
    int n = b.size();
    std::vector<double> x(n);

    c[0] /= b[0];
    d[0] /= b[0];

    for (int I = 1; I < n; I++) {
        double m = b[I] - a[I-1] * c[I-1];
        if (I < n - 1) c[I] /= m;
        d[I] = (d[I] - a[I-1] * d[I-1]) / m;
    }

    x[n-1] = d[n-1];
    for (int I = n - 2; I >= 0; I--) {
        x[I] = d[I] - c[I] * x[I+1];
    }

    return x;
}

// ****************** SHAPE FUNCTION *******************
double Solver::Rphi_1D(double ksi, int i)
{
    //return m==0 ? 0.5*(1.0-ksi) : 0.5*(ksi+1);
    double phi=1.0;
    double elementStep=2.0/m_shapeFunctionDeg; //   |----|----|----| p+1 nodes, dim(phi)=p
    double ksi_i = -1.0 + elementStep*i;
    for (int j=0; j<m_shapeFunctionDeg+1; j++)
    {
        if (i==j) continue;

        double ksi_j=-1.0+elementStep*j;
        phi*=(ksi-ksi_j)/(ksi_i-ksi_j);

    }

    return phi;
}
double Solver::Rphi_1D_deriv(double ksi, int i)
{
    //return m==0 ? -0.5 : 0.5;
    double dphi=0.0;
    double elementStep=2.0/m_shapeFunctionDeg; //   |----|----|----| p+1 nodes, dim(phi)=p
    double ksi_i = -1.0 + elementStep*i;
    for (int k=0; k<m_shapeFunctionDeg+1; k++)
    {
        if (k==i) continue;

        double ksi_k=-1.0+elementStep*k;
        double term = 1.0 / (ksi_i - ksi_k);
        for (int j=0; j<m_shapeFunctionDeg+1; j++)
        {
            if(j==k||j==i) continue;

            double ksi_j=-1.0+elementStep*j;
            term*= (ksi-ksi_j)/(ksi_i-ksi_j);
        }
        dphi += term;
    }

    return dphi;
}

// ****************************************************

// ****************** LOCAL ELEMENT *******************

void Solver::local_1D_StifnessMatrix(std::vector<double> &S, int m)
{
    double Jm=(m_grid.get_NodPosition(m+1)-m_grid.get_NodPosition(m))/2.0;
    int p = m_shapeFunctionDeg+1;
    for (int j=0; j<p; j++)
    {
        for (int i=0; i<p; i++)
        {
            int l = i+p*j;

            for (int k=0; k<m_integrationOrder; k++)
            {
                double ksi = m_local_nodes[k]; 
                double x_ksi = 0.5*((m_grid.get_NodPosition(m+1)-m_grid.get_NodPosition(m))*ksi + (m_grid.get_NodPosition(m+1)+m_grid.get_NodPosition(m)));
                S[l]+=m_local_weights[k]*(-Rphi_1D_deriv(ksi, i)*Rphi_1D_deriv(ksi,j)/Jm + B(x_ksi)*Rphi_1D_deriv(ksi,i)*Rphi_1D(ksi,j)+C(x_ksi)*Rphi_1D(ksi,i)*Rphi_1D(ksi,j)*Jm);
            }

        }
    }
}


void Solver::local_1D_LoadVector(std::vector<double> &F, int m)
{
    double Jm=(m_grid.get_NodPosition(m+1)-m_grid.get_NodPosition(m))/2.0;
    int p = m_shapeFunctionDeg+1;
    for (int j=0; j<p; j++)
    {
        for (int k=0; k<m_integrationOrder; k++)
        {
            double ksi = m_local_nodes[k]; 
            double x_ksi = 0.5*((m_grid.get_NodPosition(m+1)-m_grid.get_NodPosition(m))*ksi + (m_grid.get_NodPosition(m+1)+m_grid.get_NodPosition(m)));
            F[j]+=m_local_weights[k]*Jm*D(x_ksi)*Rphi_1D(ksi,j);
        }
    }
}


void Solver::local_1D_MassMatrix(std::vector<double> &M, int m)
{
    double Jm=(m_grid.get_NodPosition(m+1)-m_grid.get_NodPosition(m))/2.0;
    int p = m_shapeFunctionDeg+1;
    for (int j=0; j<p; j++)
    {
        for (int i=0; i<p; i++)
        {
            int l = i+p*j;

            for (int k=0; k<m_integrationOrder; k++)
            {
                double ksi = m_local_nodes[k]; 
                double x_ksi = 0.5*((m_grid.get_NodPosition(m+1)-m_grid.get_NodPosition(m))*ksi + (m_grid.get_NodPosition(m+1)+m_grid.get_NodPosition(m)));
                M[l]+=m_local_weights[k]*(Rphi_1D(ksi,i)*Rphi_1D(ksi,j)*Jm);
            }

        }
    }
}
// ****************************************************

// ****************** BOUNDARY CONDITIONS *******************

void Solver::boundaryConditions(std::string work_type)
{
    int p = m_shapeFunctionDeg+1;
    int N = (p-1)*m_grid.get_elementNumber()+1;

    std::cout<<"Including boundary conditions...\n";
    if (work_type=="stationary")
    {
        for (const auto &bc : m_grid.get_BC())
        {
            if (bc.bct==1)
            {
                int m = bc.m;
                int I_glob = m*(p-1);
                for (int j=0; j<N; j++)
                {
                    int L_glob = I_glob + j*N;

                    m_S1D[L_glob]=0.0;

                    if(I_glob==j)
                        m_S1D[L_glob]=1.0;
                }

                m_F1D[I_glob]=bc.bc_value;
            }
        }
    }

    if (work_type=="eigen")
    {
        for (const auto &bc : m_grid.get_BC())
        {
            if (bc.bct==1)
            {
                int m = bc.m;
                int I_glob = m*(p-1);
                for (int j=0; j<N; j++)
                {
                    int L_glob = I_glob + j*N;

                    m_S1D[L_glob]=0.0;
                    m_M1D[L_glob]=0.0;

                    if(I_glob==j)
                    {
                        m_S1D[L_glob]=1.0;
                        m_M1D[L_glob]=1.0;
                    }
                }    


            }
        }
    }
    std::cout<<"Including boundary conditions DONE\n";

}


// ****************************************************

// ****************** ASSEMBLERS *******************

void Solver::Matrix_assembler(std::string work_type)
{
    int p = m_shapeFunctionDeg+1;
    std::vector<double> S_loc(p*p, 0.0);
    int N = (p-1)*m_grid.get_elementNumber()+1; 

    m_S1D.resize(N*N, 0.0);

    std::cout<<"Stifness assemble process start...\n";

    for (int m=0; m<m_grid.get_elementNumber(); m++)
    {
        S_loc.assign(S_loc.size(), 0.0);

        local_1D_StifnessMatrix(S_loc, m);

        for (int j=0; j<p; j++)
        {
            int J_glob = m*(p-1)+j;
            for(int i=0; i<p; i++)
            {
                int I_glob=m*(p-1)+i;

                int L_glob=I_glob+N*J_glob;
                int l = i+p*j;
                m_S1D[L_glob]+=S_loc[l];
            }
        }
    }

    std::cout<<"Stfiness assemble process DONE\n";


    if (work_type=="eigen")
    {
        std::cout<<"Mass assemble process start...\n";

        std::vector<double> M_loc(p*p, 0.0);
        m_M1D.resize(N*N, 0.0);
        for (int m=0; m<m_grid.get_elementNumber(); m++)
        {
            M_loc.assign(M_loc.size(), 0.0);

            local_1D_MassMatrix(M_loc, m);

            for (int j=0; j<p; j++)
            {
                int J_glob = m*(p-1)+j;
                for(int i=0; i<p; i++)
                {
                    int I_glob=m*(p-1)+i;

                    int L_glob=I_glob+N*J_glob;
                    int l = i+p*j;
                    m_M1D[L_glob]+=M_loc[l];
                }
            }
        }

    std::cout<<"Mass assemble process DONE\n";
    }

}

void Solver::Vector_assembler()
{
    int p = m_shapeFunctionDeg+1;
    std::vector<double> F_loc(p, 0.0);
    int N = (p-1)*m_grid.get_elementNumber()+1; 
    m_F1D.resize(N, 0.0);

    std::cout<<"Load assemble process start...\n";

    for (int m=0; m<m_grid.get_elementNumber(); m++)
    {
        F_loc.assign(F_loc.size(), 0.0);

        local_1D_LoadVector(F_loc, m);

        for (int j=0; j<p; j++)
        {
            int J_glob = m*(p-1)+j;
            m_F1D[J_glob]+=F_loc[j];
        }
    }

    std::cout<<"Load assemble process DONE\n";

}



// ****************************************************

// ****************** SOLVER STATIONARY 1 DIM *******************

void Solver::stationary_1D()
{
    std::cout<<"Solving System of Eq...\n";
    
    int N = m_shapeFunctionDeg*m_grid.get_elementNumber()+1; 
    m_Q.resize(N);
    Eigen::Map<Eigen::MatrixXd> S_e(m_S1D.data(), N, N);
    Eigen::Map<Eigen::VectorXd> F_e(m_F1D.data(), N);
    Eigen::Map<Eigen::VectorXd> Q_e(m_Q.data(), N);

    Q_e = S_e.partialPivLu().solve(F_e);

    std::cout<<"Solving System of Eq DONE\n";
}

void Solver::Eigen_1D()
{
    std::cout<<"Eigen solver start..\n"; 
    int N = m_shapeFunctionDeg*m_grid.get_elementNumber()+1; 
    m_Q.resize(N);
    Eigen::Map<Eigen::MatrixXd> S_e(m_S1D.data(), N, N);
    Eigen::Map<Eigen::MatrixXd> M_e(m_M1D.data(), N, N);

    Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> solver(S_e, M_e);

    Eigen::VectorXd realEigenvalues = solver.eigenvalues().real();
    Eigen::MatrixXd realEigenvectors = solver.eigenvectors().real();

    std::vector<std::pair<double, int>> eigenPairs;
    for (int i = 0; i < N; ++i) {
        eigenPairs.push_back({realEigenvalues(i), i});
    }

    std::sort(eigenPairs.begin(), eigenPairs.end(), [](const auto& a, const auto& b) {
        return a.first > b.first; 
    });

    m_Q.resize(N);
    m_eigenValues1D.resize(N);
    m_eigenVectors1D.resize(N * N);

    for (int mode = 0; mode < N; ++mode) {
        m_eigenValues1D[mode] = eigenPairs[mode].first;
        
        int orig_mode = eigenPairs[mode].second;
        
        for (int node = 0; node < N; ++node) {

            m_eigenVectors1D[node + mode * N] = realEigenvectors(node, orig_mode);
        }
    }

    std::cout<<"Eigen solver DONE\n";
}

//==============================================

//================== SOLUTION ==================

double Solver::U_1D(double x)
{
    int N = m_grid.get_elementNumber();
    int p_deg = m_shapeFunctionDeg;

    for (int m = 0; m < N; m++)
    {
        double x_left = m_grid.get_NodPosition(m);
        double x_right = m_grid.get_NodPosition(m+1);

        if (x >= x_left && x <= x_right)
        {

            double ksi = (2.0 * (x - x_left) / (x_right - x_left)) - 1.0;
            
            double u_val = 0.0;
            for (int i = 0; i <= p_deg; i++)
            {
                int glob_idx = m * p_deg + i; 
                u_val += m_Q[glob_idx] * Rphi_1D(ksi, i);
            }
            
            return u_val;
        }
    }
    return 0.0;
}

double Solver::U_1D(double x, int mode)
{
    int N = m_grid.get_elementNumber();
    int p_deg = m_shapeFunctionDeg;

    for (int m = 0; m < N; m++)
    {
        double x_left = m_grid.get_NodPosition(m);
        double x_right = m_grid.get_NodPosition(m+1);

        if (x >= x_left && x <= x_right)
        {
            double ksi = (2.0 * (x - x_left) / (x_right - x_left)) - 1.0;
            
            double u_val = 0.0;
            for (int i = 0; i <= p_deg; i++)
            {
                int glob_idx = m * p_deg + i;

                int l = glob_idx + mode*(p_deg*m_grid.get_elementNumber()+1); 
                u_val += m_eigenVectors1D[l] * Rphi_1D(ksi, i);
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
        std::cout<<"Saving data...\n";
        auto [x_start, x_end] = m_grid.get_RodInfo();
        double x=x_start;

        std::ofstream file(outdir+"/u.dat");

        while (x<=x_end)
        {
            file<<x<<" "<<U_1D(x)<<"\n";
            x+=m_dxSave;
        }
    
        file<<x_end<<" "<<U_1D(x_end)<<"\n";

        file.close();
        std::cout<<"Saving data DONE\n";
    }

    else if (work_type=="eigen")
    {
        
        std::ofstream file_eigenvalue(outdir+"/eigenvalue.dat");
        
        for (int mode=0; mode<m_shapeFunctionDeg*m_grid.get_elementNumber()+1; mode++)
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









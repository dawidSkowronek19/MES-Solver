#include "MathSolver.h"

//==================== EIGEN SECTION ==============

std::pair<Eigen::VectorXd, Eigen::MatrixXd> MathSolver::EigenSymetric_real(std::vector<double> &S, std::vector<double> &M)
{
    std::cout<<"\tSOLVER METHOD: ";
    std::cout<<"Eigen, Self-Adjont decomposition\n";
    int N = sqrt(S.size());
    
    Eigen::Map<Eigen::MatrixXd> S_e(S.data(), N, N);
    Eigen::Map<Eigen::MatrixXd> M_e(M.data(), N, N);

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(S_e, M_e);

    //std::cout<<solver.info()<<"\n";

    Eigen::VectorXd Eigenvalues = solver.eigenvalues();
    Eigen::MatrixXd Eigenvectors = solver.eigenvectors();


return {Eigenvalues, Eigenvectors};

}


std::pair<Eigen::VectorXd, Eigen::MatrixXd> MathSolver::EigenGeneral_real(std::vector<double> &S, std::vector<double> &M)
{
    std::cout<<"\tSOLVER METHOD: ";
    std::cout<<"Eigen, General solver\n";
    int N = sqrt(S.size());
    
    Eigen::Map<Eigen::MatrixXd> S_e(S.data(), N, N);
    Eigen::Map<Eigen::MatrixXd> M_e(M.data(), N, N);

    //std::cout<<"Mapping 2 done\n";

    Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> solver(S_e, M_e);

    //std::cout<<solver.info()<<"\n";

    double max_imaginary = solver.eigenvalues().imag().cwiseAbs().maxCoeff();
    if (max_imaginary > 1e-10) 
    {
        throw std::runtime_error("FATAL ERROR: Complex Eigen Values\n");
    }

    Eigen::VectorXd realEigenvalues = solver.eigenvalues().real();
    Eigen::MatrixXd realEigenvectors = solver.eigenvectors().real();

    std::vector<std::pair<double, int>> eigenPairs;
    for (int i = 0; i < N; ++i) {
        eigenPairs.push_back({realEigenvalues(i), i});
    }
    
    std::sort(eigenPairs.begin(), eigenPairs.end(), [](const auto& a, const auto& b) {
        return a.first < b.first; 
    });
    
    Eigen::MatrixXd sorted_realEigenVectors(N,N);
    Eigen::VectorXd sorted_realEigenValues(N);

    for (int mode = 0; mode < N; ++mode) {
        sorted_realEigenValues[mode] = eigenPairs[mode].first;
        
        int orig_mode = eigenPairs[mode].second;
        
        for (int node = 0; node < N; ++node) {

            sorted_realEigenVectors(node, mode) = realEigenvectors(node, orig_mode);
        }
    }
    return {sorted_realEigenValues, sorted_realEigenVectors};
}


// ================ System of Equations Section ==============

Eigen::VectorXd MathSolver::StEqSymetric_real(std::vector<double> &S, std::vector<double> &F)
{
    std::cout<<"\tSOLVER METHOD: ";
    std::cout<<"System of Equations, Bunch-Kaufman decomposition\n";
    int N=F.size();
    Eigen::Map<Eigen::MatrixXd> S_e(S.data(), N, N);
    Eigen::Map<Eigen::VectorXd> F_e(F.data(), N);

    
    Eigen::LDLT<Eigen::MatrixXd> solver;
    solver.compute(S_e);
    
    //std::cout<<solver.info()<<"\n";
    return solver.solve(F_e);
}

Eigen::VectorXd MathSolver::StEqGeneral_real(std::vector<double> &S, std::vector<double> &F)
{
    std::cout<<"\tSOLVER METHOD: ";
    std::cout<<"System of Equations, LU decomposition\n";
    int N=F.size();
    Eigen::Map<Eigen::MatrixXd> S_e(S.data(), N, N);
    Eigen::Map<Eigen::VectorXd> F_e(F.data(), N);

    return S_e.partialPivLu().solve(F_e);
}


#include "SolverPhysics.h"

//==================== PHYSICS PART ==============================

// *****************************************************************
Physics::Physics(MathSolver &math) : m_math(math)
{
    std::cout<<std::string(60, '=')<<"\n";
    std::cout<<"\t\t\tPHYSICS SECTION\n\n";
    std::cout<<"EQUATION NAME: ";

}


std::pair<Eigen::VectorXd, Eigen::MatrixXd> Physics::EigenSolver(std::vector<double> &S, std::vector<double> &M)
{
    std::cout<<"NO EIGEN PROBLEM FOR THIS EQUATION\n";

    return {};
}

double Physics::B(double x, double u) {return 0.0;}
double Physics::C(double x, double u) {return 0.0;}
double Physics::D(double x, double u) {return 0.0;}




// *****************************************************************

//Poisson

Poisson::Poisson(MathSolver &math) : Physics(math)
{
    std::cout<<"Poisson\n";
}

double Poisson::D(double x, double u) {return -1.0;}


Eigen::VectorXd Poisson::solver(std::vector<double> &S, std::vector<double> &F)
{
    
    return m_math.StEqSymetric_real(S, F);
}



//General PDE

GeneralPDE::GeneralPDE(MathSolver &math) : Physics(math)
{
    std::cout<<"General form PDE\n";
}

double GeneralPDE::B(double x, double u) {return 0.0;}
double GeneralPDE::C(double x, double u) {return 0.0;}
double GeneralPDE::D(double x, double u) {return 0.0;}


Eigen::VectorXd GeneralPDE::solver(std::vector<double> &S, std::vector<double> &F)
{
    return m_math.StEqGeneral_real(S,F);
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> GeneralPDE::EigenSolver(std::vector<double> &S, std::vector<double> &M)
{
    return m_math.EigenGeneral_real(S, M);
}

//General symetric PDE

GeneralSymetricPDE::GeneralSymetricPDE(MathSolver &math) : Physics(math)
{
    std::cout<<"General symetric form PDE\n";
}

double GeneralSymetricPDE::C(double x, double u) {return 0.0;}
double GeneralSymetricPDE::D(double x, double u) {return 2.0*u*u*u;}



Eigen::VectorXd GeneralSymetricPDE::solver(std::vector<double> &S, std::vector<double> &F)
{
    return m_math.StEqSymetric_real(S,F);
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> GeneralSymetricPDE::EigenSolver(std::vector<double> &S, std::vector<double> &M)
{
    return m_math.EigenSymetric_real(S,M);
}


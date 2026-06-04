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

double Physics::B(double x) {return 0.0;}
double Physics::C(double x) {return 0.0;}
double Physics::D(double x) {return 0.0;}

double Physics::dB_du(double x, double u) {return 0.0;}
double Physics::dC_du(double x, double u) {return 0.0;}
double Physics::dD_du(double x, double u) {return 0.0;}




// *****************************************************************

//Poisson

Poisson::Poisson(MathSolver &math) : Physics(math)
{
    std::cout<<"Poisson\n";
}

double Poisson::D(double x) {return -1.0;}


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

double GeneralPDE::B(double x) {return 0.0;}
double GeneralPDE::C(double x) {return 0.0;}
double GeneralPDE::D(double x) {return 0.0;}

double GeneralPDE::dB_du(double x, double u) {return 0.0;}
double GeneralPDE::dC_du(double x, double u) {return 0.0;}
double GeneralPDE::dD_du(double x, double u) {return 0.0;}


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
double GeneralSymetricPDE::D(double x, double u) {return 2*u*u*u;}

double GeneralSymetricPDE::dC_du(double x, double u) {return 0.0;}
double GeneralSymetricPDE::dD_du(double x, double u) {return 6*u*u;}

double GeneralSymetricPDE::C(double x) {return 0.0;}
double GeneralSymetricPDE::D(double x) {return 0.0;}



Eigen::VectorXd GeneralSymetricPDE::solver(std::vector<double> &S, std::vector<double> &F)
{
    return m_math.StEqSymetric_real(S,F);
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> GeneralSymetricPDE::EigenSolver(std::vector<double> &S, std::vector<double> &M)
{
    return m_math.EigenSymetric_real(S,M);
}


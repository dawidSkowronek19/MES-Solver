#include "SolverPhysics.h"

//==================== PHYSICS PART ==============================

// *****************************************************************
Physics::Physics(MathSolver &math, ConfigParameters &config) : m_math(math), m_config(config)
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
double Physics::E(double x) {return 0.0;}
double Physics::F(double x) {return 0.0;}


double Physics::dB_du(double x, double u) {return 0.0;}
double Physics::dC_du(double x, double u) {return 0.0;}
double Physics::dD_du(double x, double u) {return 0.0;}


void Physics::connectTime(const double *t) 
{
    m_t=t;
}

double Physics::getCurrentTime(){return (m_t!=nullptr) ? *m_t : 0.0; }

// *****************************************************************

//Poisson

Poisson::Poisson(MathSolver &math, ConfigParameters &config) : Physics(math, config)
{
    std::cout<<"Poisson\n";
}

double Poisson::D(double x) 
{
    return m_config.Functions[2](x,0.0, 0.0);
}


Eigen::VectorXd Poisson::solver(std::vector<double> &S, std::vector<double> &F)
{
    
    return m_math.StEqSymetric_real(S, F);
}



//General PDE

GeneralPDE::GeneralPDE(MathSolver &math, ConfigParameters &config) : Physics(math, config)
{
    std::cout<<"General form PDE\n";
}

double GeneralPDE::B(double x, double u) 
{
    double t=getCurrentTime();
    return m_config.Functions[0](x,u,t);
}
double GeneralPDE::C(double x, double u) 
{
    double t=getCurrentTime();
    return m_config.Functions[1](x,u,t);
}
double GeneralPDE::D(double x, double u) 
{
    double t=getCurrentTime();
    return m_config.Functions[2](x,u,t);
}

double GeneralPDE::B(double x) 
{
    double t=getCurrentTime();
    return m_config.Functions[0](x,0.0,t);
}
double GeneralPDE::C(double x) 
{
    double t=getCurrentTime();
    return m_config.Functions[1](x,0.0,t);
}
double GeneralPDE::D(double x) 
{
    double t=getCurrentTime();
    return m_config.Functions[2](x,0.0,t);
}
double GeneralPDE::E(double x) 
{
    double t=getCurrentTime();
    return m_config.Functions[3](x,0.0,t);
}
double GeneralPDE::F(double x) 
{
    double t=getCurrentTime();
    return m_config.Functions[4](x,0.0,t); //finish
}

double GeneralPDE::dB_du(double x, double u) 
{
    double t=getCurrentTime();
    return (B(x,u+du)-B(x,u-du))/(2.0*du);
}
double GeneralPDE::dC_du(double x, double u) 
{
    double t=getCurrentTime();
    return (C(x,u+du)-C(x,u-du))/(2.0*du);;
}
double GeneralPDE::dD_du(double x, double u) 
{
    double t=getCurrentTime();
    return (D(x,u+du)-D(x,u-du))/(2.0*du);;
}


Eigen::VectorXd GeneralPDE::solver(std::vector<double> &S, std::vector<double> &F)
{
    return m_math.StEqGeneral_real(S,F);
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> GeneralPDE::EigenSolver(std::vector<double> &S, std::vector<double> &M)
{
    return m_math.EigenGeneral_real(S, M);
}

//General symetric PDE

GeneralSymetricPDE::GeneralSymetricPDE(MathSolver &math, ConfigParameters &config) : Physics(math, config)
{
    std::cout<<"General symetric form PDE\n";
}

double GeneralSymetricPDE::C(double x, double u) 
{
    double t=getCurrentTime();
    return m_config.Functions[1](x,u,t);
}
double GeneralSymetricPDE::D(double x, double u) 
{
    double t=getCurrentTime();
    return m_config.Functions[2](x,u,t);
}

double GeneralSymetricPDE::dC_du(double x, double u) 
{
    double t=getCurrentTime();
    return (C(x,u+du)-C(x,u-du))/(2.0*du);
}
double GeneralSymetricPDE::dD_du(double x, double u) 
{
    double t=getCurrentTime();
    return (D(x,u+du)-D(x,u-du))/(2.0*du);
}

double GeneralSymetricPDE::C(double x) 
{
    double t=getCurrentTime();
    return m_config.Functions[1](x,0.0,t);
}
double GeneralSymetricPDE::D(double x) 
{
    double t=getCurrentTime();
    return m_config.Functions[2](x,0.0,t);
}
double GeneralSymetricPDE::E(double x) 
{
    double t=getCurrentTime();
    return m_config.Functions[3](x,0.0,t);
}
double GeneralSymetricPDE::F(double x) 
{
    double t=getCurrentTime();
    return m_config.Functions[4](x,0.0,t);
}



Eigen::VectorXd GeneralSymetricPDE::solver(std::vector<double> &S, std::vector<double> &F)
{
    return m_math.StEqSymetric_real(S,F);
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> GeneralSymetricPDE::EigenSolver(std::vector<double> &S, std::vector<double> &M)
{
    return m_math.EigenSymetric_real(S,M);
}


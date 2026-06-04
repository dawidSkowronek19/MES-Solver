#include "ShapeFunction.h"

// =================== Shape Function Part =======================

ShapeFunction::ShapeFunction(int shapeFunctionDeg) : m_shapeFunctionDeg(shapeFunctionDeg)
{
    //std::cout<<"ShapeFunctionCreated\n";
}

ShapeFunction::~ShapeFunction() {}

double ShapeFunction::Rphi_1D(double ksi, int i)
{
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

double ShapeFunction::Rphi_1D_deriv(double ksi, int i)
{
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

int ShapeFunction::get_deg() {return m_shapeFunctionDeg;}
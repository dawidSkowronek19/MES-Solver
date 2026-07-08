#include "shapeFunction.hpp"

ShapeFunction::ShapeFunction(const ElementPointPositions &acc_element, int p) : m_accElement(acc_element), m_p(p)
{
    for (int k=0; k<=m_p; k++)
    {
        for (int j=0; j<=m_p-k; j++)
        {
            int i = m_p - j - k;
            m_IJK.push_back({i,j,k});
        }
    }
}


//============================== SHAPE FUNCTION INIT ===================================

std::tuple<double, double, double> ShapeFunction::lambda_gen(const double ksi, const double eta) const
{
    return {1.0-ksi-eta, ksi, eta};
}

double ShapeFunction::Silvester(const int i, const double lambda) const
{
    double value=1.0;
    for (int idx=0; idx<i; idx++)
    {

        double T= m_p*lambda - idx; 
        double B= i - idx;

        value*=(T/B);
    }

    return value;
}

double ShapeFunction::divSilvester(const int i, const double lambda) const
{
    double div=0.0;
    for (int sum_idx=0; sum_idx<i; sum_idx++)
    {
        double mult=1.0;
        for (int mult_idx=0; mult_idx<i; mult_idx++)
        {
            if (mult_idx==sum_idx)  continue;

            mult*=(m_p*lambda-mult_idx)/(i-mult_idx);
        }

        div+=(static_cast<double>(m_p) / (i - sum_idx)) * mult;
    }

    return div;
}



double ShapeFunction::phi(double ksi, double eta, int idx) const
{
    auto [lambda1, lambda2, lambda3] = lambda_gen(ksi, eta);
    auto [I,J,K] = m_IJK[idx];

    return Silvester(I,lambda1)*Silvester(J,lambda2)*Silvester(K, lambda3);
}

std::pair<double, double> ShapeFunction::div_phi(double ksi, double eta, int idx) const
{
    double dphi_dksi;
    double dphi_deta;

    auto [lambda1, lambda2, lambda3] = lambda_gen(ksi, eta);
    auto [I,J,K] = m_IJK[idx];

    double S_I = Silvester(I, lambda1);
    double S_J = Silvester(J, lambda2);
    double S_K = Silvester(K, lambda3);

    double dlambda1_dksi = -1.0;
    double dlambda2_dksi = 1.0;
    double dlambda1_deta = -1.0;
    double dlambda3_deta=1.0; 

    dphi_dksi = divSilvester(I,lambda1)*dlambda1_dksi*S_J*S_K + S_I*divSilvester(J,lambda2)*dlambda2_dksi*S_K;
    dphi_deta = divSilvester(I,lambda1)*dlambda1_deta*S_J*S_K + S_I*S_J*divSilvester(K,lambda3)*dlambda3_deta;

    return {dphi_dksi, dphi_deta};

}


void ShapeFunction::get_cached_values(const std::vector<std::pair<double, double>>&integration_points)
{
    int shfunc_nb = m_IJK.size();
    
    for (int shfunc_idx=0; shfunc_idx<shfunc_nb; shfunc_idx++)
    {
        for (int point_idx=0; point_idx<integration_points.size(); point_idx++)
        {
            const auto [ksi, eta] = integration_points[point_idx];
            const auto [dphi_dksi, dphi_deta] = div_phi(ksi, eta, shfunc_idx);
            m_phi(point_idx, shfunc_idx) = phi(ksi, eta, shfunc_idx);
            m_dphi_dksi(point_idx, shfunc_idx) = dphi_dksi;
            m_dphi_deta(point_idx, shfunc_idx) = dphi_deta;
        }
    }
}

double ShapeFunction::get_phi(const int point_idx, const int sh_nb) const {return m_phi(point_idx, sh_nb);}
double ShapeFunction::get_dphi_dksi(const int point_idx, const int sh_nb) const {return m_dphi_dksi(point_idx, sh_nb);}
double ShapeFunction::get_dphi_deta(const int point_idx, const int sh_nb) const {return m_dphi_deta(point_idx, sh_nb);}
int ShapeFunction::get_p() const {return m_p;}
//===================================================================================



// ==================================================================================


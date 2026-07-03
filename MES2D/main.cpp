#include "./mesh/Grid.hpp"
#include "./solver/shapeFunction.hpp"

#include <fstream>

int main()
{
    GeometryParameters geoparam{
        100,
        2*sin(M_PI/100),
        2*0.5*sin(M_PI/100),
        100
    };
    int p = 2;
    Grid2D testgrid(geoparam);
    testgrid.triangular();
    testgrid.create_neighboursList();
    testgrid.relaxGrid();
    testgrid.savePointsList();
    testgrid.saveTrianglesPoints();

    ShapeFunction sh_func(testgrid, p);

    int N_ksi = 100;
    int N_eta=100;
    double dksi = 1.0/(N_ksi-1);
    double deta = 1.0/(N_eta-1);
    int N = (p+1)*(p+2)/2;
    for (int p_idx=0; p_idx<N; p_idx++)
    {
        std::ofstream file_phi("./outdir/phi_"+std::to_string(p_idx)+".dat");
        std::ofstream file_dphi("./outdir/dphi_"+std::to_string(p_idx)+".dat");
        for (int ksi_idx=0; ksi_idx<N_ksi; ksi_idx++)
        {
            double ksi = ksi_idx*dksi;
            for (int eta_idx=0; eta_idx<N_eta-ksi_idx; eta_idx++)
            {
                double eta = eta_idx*deta;
                auto [dphi_dksi, dphi_deta] = sh_func.div_phi(ksi, eta, p_idx);

                file_phi<<ksi<<" "<<eta<<" "<<sh_func.phi(ksi, eta, p_idx)<<"\n";
                file_dphi<<ksi<<" "<<eta<<" "<<dphi_dksi<<" "<<dphi_deta<<"\n";
            }
        }
        file_phi.close();
        file_dphi.close();
    }


}   

#include "Grid.h"

//========== Constructor/Destructor ==========

Grid_1D::Grid_1D(double x_start, double x_end, double growFactor, double A, int int_order, int numbOfElements): m_x_start(x_start), m_x_end(x_end), 
            m_growFactor(growFactor), m_A(A), m_N_int(int_order), m_numbOfElements(numbOfElements)
            {
                
                std::cout<<std::string(60, '=');
                
                std::cout<<"\n\n\n \t\t\tGRID SECTION \n\n\n";
                
                std::cout<<"\t\t"<<std::string(4, '#')<<" GRID PARAMETERS "<<std::string(4, '#')<<"\n\n";
                std::cout<<"\tx_start = "<<m_x_start<<"\n";
                std::cout<<"\tx_end = "<<m_x_end<<"\n";
                std::cout<<"\tGrow Factor = "<<m_growFactor<<"\n";
                std::cout<<"\tAttraction Constant = "<<m_A<<"\n";
                std::cout<<"\tGrid Integration Order = "<<m_N_int<<"\n";
                std::cout<<"\tNumber of elements = "<<m_numbOfElements<<"\n\n";
                std::cout<<"\t\t"<<"# GRID PARAMETERS END #\n\n";
                m_nodes.resize(m_numbOfElements+1, 0.0);
            }

Grid_1D::~Grid_1D(){}
//============================================

//========= Setting Boundary conditions ======
void Grid_1D::set_boundaryConditions(std::vector<BoundaryCondition> bc)
{
    std::cout<<"READING BOUNDARY CONDITIONS FROM INPUT FILE\n";

    //0 - no condition
    //1 - Dirichlet
    //2 - Neuman
    //3 -Robin
    //4- Born

    m_bc=bc;


    std::cout<<std::string(4, '#')<<" BOUNDARY CONDITIONS "<<std::string(4, '#')<<"\n\n";

    for (auto &bc : m_bc)
    {
        if (bc.bct==1)
        {
            std::cout<<"\t Dirichlet Condition ";
        }
        else if (bc.bct==2)
        {
            std::cout<<"\t Neumann Condition ";
        }

        std::cout<<"\tValue = "<<bc.bc_value<<"\n";
        std::cout<<"\tCondition set in element nb "<<bc.m<<"\n\n";

    }

    std::cout<<"# BOUNDARY CONDITIONS END #\n\n";
    
}
//===========================================

//============ Build Grid ==================

void Grid_1D::buildGrid()
{

    //building node need function
    std::cout<<"BUILDING GRID...";
    int N=m_bc.size();
    double ev_dx=(m_x_end-m_x_start)/m_numbOfElements;
    std::vector<double> boundCondition_pos;

    for (int m=0; m<N; m++)
    {
        boundCondition_pos.push_back(m_bc[m].m*ev_dx + m_x_start);
    }
    

    double sigma=m_growFactor;
    double A=m_A;
    auto node_need = [N, boundCondition_pos, sigma, A](double x) -> double
    {
        double omega=1.0;
        for (int i=0; i<N; i++)
        {
            omega+=A*exp(-(x-boundCondition_pos[i])*(x-boundCondition_pos[i])/(2.0*sigma*sigma));
        }

        return omega;
    };
    // computing distribuante
    //std::cout<<"Computing dist...\n";
    const gsl_integration_fixed_type * T = gsl_integration_fixed_legendre;
    
    std::vector<double> F(m_numbOfElements+1, 0.0);
    for (int m=1; m<m_numbOfElements+1; m++)
    {
        double x=m_x_start+m*ev_dx;
        int integration_order=m+m_N_int;
        gsl_integration_fixed_workspace *work = gsl_integration_fixed_alloc(T, integration_order, m_x_start, x, 0, 0);
        double *int_nods=gsl_integration_fixed_nodes(work);
        double *int_w=gsl_integration_fixed_weights(work);

        for (int x_idx=0; x_idx<integration_order; x_idx++)
        {
            F[m]+=int_w[x_idx]*node_need(int_nods[x_idx]);
        }

        gsl_integration_fixed_free(work);
    }
    
    //computing nodes
    //std::cout<<"Computing nodes positions...\n";
    double F_L=F.back();
    for (int j=0; j<m_numbOfElements+1; j++)
    {
        double targetF=j*F_L/m_numbOfElements;
        auto it = std::lower_bound(F.begin(), F.end(), targetF);

        if (it==F.begin())
        {
            m_nodes[j]=m_x_start;
        }
        else
        {
            size_t idx1 = std::distance(F.begin(), it);
            size_t idx0 = idx1 - 1;

            double F0 = F[idx0], F1 = F[idx1];
            double x0 = m_x_start + idx0 * ev_dx;
            double x1 = m_x_start + idx1 * ev_dx;

            double t = (targetF - F0) / (F1 - F0);
            m_nodes[j] = x0 + t * (x1 - x0);
        }
    }

    std::cout<<"DONE\n";


}
//==========================================

//================ SAVING GRID =============

void Grid_1D::saveGrid(std::string outdir)
{
    std::ofstream file(outdir+"/mesh.dat");

    file<<"node number\tposition\n";

    for (int m=0; m<m_numbOfElements+1; m++)
    {
        file<<m<<" "<<m_nodes[m]<<"\n";
    }

    file.close();

    std::cout<<"\n\n\n \t\t\tGRID SECTION DONE\n";
}

//===========================================

//========= Getters ================
double Grid_1D::get_NodPosition(int m) const {return m_nodes[m];}

double Grid_1D::get_elementNumber() const {return m_numbOfElements;}

std::pair<double, double> Grid_1D::get_RodInfo() const {return {m_x_start, m_x_end};}

std::vector<BoundaryCondition> Grid_1D::get_BC() const {return m_bc;}
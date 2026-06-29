#ifndef SHAPE_FUNC
#define SHAPE_FUNC
#include<iostream>
#include "Parser.h"
class ShapeFunction{

    public:
        ShapeFunction(ConfigParameters &config);
        ~ShapeFunction();
        double Rphi_1D(double ksi, int i);
        double Rphi_1D_deriv(double ksi, int i);

        int get_deg();
    private:
        ConfigParameters m_config;


};

#endif
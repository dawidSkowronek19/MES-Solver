#ifndef SHAPE_FUNC
#define SHAPE_FUNC
#include<iostream>
class ShapeFunction{

    public:
        ShapeFunction(int shapeFunctionDeg);
        ~ShapeFunction();
        double Rphi_1D(double ksi, int i);
        double Rphi_1D_deriv(double ksi, int i);

        int get_deg();
    private:
        int m_shapeFunctionDeg;


};

#endif
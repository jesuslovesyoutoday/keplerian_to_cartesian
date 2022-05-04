#include "kepl_to_cart.h"
#include <iostream>

int main() {
    
    double a = 57909226541.52;
    double e = 0.20563593;
    double omega = 77.45779628 * 3.14/180; //1.35
    double big_omega = 48.33076593 * 3.14 / 180; // 0.84
    double i = 7.00497902 * 3.14/180; // 0.1221
    double m0 = 252.25032350 * 3.14/180; // 4.76
    //double mu = 3.986 * 10e14 ;
    //double mu = 1.327124 * 10e20;
    double mu = 1.98e+30;
    
    BodyPosition obj0 = BodyPosition(a, e, omega, big_omega, i, m0, mu);
    obj0.forwardConvertation();
    double* r  = obj0.getPosition();
    double* dr = obj0.getVelocity();
    /*for (int i = 0; i < 3; i++)
    {
        std::cout << "v  " << r[i] << std::endl;
        std::cout << "dv " << dr[i] << std::endl;   
    }*/
    
    double  r0[3] = {5.60950066e+10, 1.84707796e+10, -3.63983421e+09};
    double dr0[3] = {-2.30603729e+04, 4.04106010e+04, 5.41772383e+03};
    
    BodyPosition obj = BodyPosition(r0, dr0, mu);
    obj.backConvertation();
    double a_ = obj.getA();
    double e_ = obj.getE();
    double omega_ = obj.getOmega();
    double big_omega_ = obj.getBigOmega();
    double i_ = obj.getI();
    double m0_ = obj.getM0();
    std::cout << "a  " << a_ << std::endl;
    std::cout << "e  " << e_ << std::endl;
    std::cout << "omega  " << omega_ << std::endl;
    std::cout << "big_omega  " << big_omega_ << std::endl;
    std::cout << "i  " << i_ << std::endl;
    std::cout << "m0  " << m0_ << std::endl;

    return 0;
}

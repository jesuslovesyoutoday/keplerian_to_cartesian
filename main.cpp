#include "kepl_to_cart.h"
#include <iostream>

//TODO: - main
//      - data from noaa
//      - readme
//      - tests (with astropy)

int main() {
    
    double a = 10e5;
    //double a = 149597870700.0;
    double e = 0.0006703;
    //double e = 0;
    double omega = 130.5360 * 3.14 / 180;
    //double omega = 0;
    //double big_omega = 0;
    double big_omega = 247.4627 * 3.14 / 180;
    double i = 51.6416 * 3.14;
    //double i = 0;
    double m0 = 1;
    //double m0 = 0;
    //double mu = 3.986 * 10e14 ;
    double mu = 1.327124 * 10e20;

    BodyPosition obj = BodyPosition(a, e, omega, big_omega, i, m0, mu);
    obj.calculate();
    double* r  = obj.getPosition();
    double* dr = obj.getVelocity();

    for (int i = 0; i < 3; i++)
    {
        std::cout << "v  " << r[i] << std::endl;
        std::cout << "dv " << dr[i] << std::endl;   
    }

    return 0;
}

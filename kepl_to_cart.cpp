#include "kepl_to_cart.h"
#include <iostream>
#include <cmath>

BodyPosition::BodyPosition(double a, double e, double omega, 
                           double big_omega, double i, double m0, double mu) 
{
    this->a         = a;
    this->e         = e;
    this->omega     = omega;
    this->big_omega = big_omega;
    this->i         = i;
    this->m0        = m0;
    this->mu        = mu;
}

void BodyPosition::calculate() {

    double x  = sqrt(1 - this->e) * cos(this->m0/2);
    double y  = sqrt(1 + this->e) * sin(this->m0/2);
    double nu = 2 * atan2(y, x); 
    std::cout << "nu" << " " << nu << std::endl;

    double rc = this->a * (1 - this->e * cos(this->m0));
    std::cout << "rc " << rc << std::endl; 
    
    double O[3]     = {rc * cos(nu), rc * sin(nu), 0};
    std::cout << "o" << O[0] << O[1] << O[2] << std::endl;
    double param    = (sqrt(nu * this->a))/rc;
    double dO_x     = - sin(this->m0);
    double dO_y     = sqrt(1 - pow(this->e, 2)) * cos(this->m0);
    double dO_z     = 0;
    double dO[3]    = {param * dO_x, param * dO_y, param * dO_z};
    std::cout << "do" << dO[0] << dO[1] << dO[2] << std::endl;


    double r_x = (O[0] * (cos(this->omega) * cos(this->big_omega)
                       - sin(this->omega) * sin(this->big_omega) 
                       * cos(this->i))
                - O[1] * (sin(this->omega) * cos(this->big_omega)
                       + cos(this->omega) * cos(this->i) 
                       * sin(this->big_omega)));
    double r_y = (O[0] * (cos(this->omega) * sin(this->big_omega)
                       + sin(this->omega) * cos(this->big_omega) 
                       * cos(this->i))
                + O[1] * (cos(this->omega) * cos(this->big_omega)
                       * cos(this->i) - sin(this->big_omega) 
                       * sin(this->omega)));
    double r_z = (O[0] * sin(this->omega) * sin(this->i)
                + O[1] * cos(this->omega) * sin(this->i));
    this->r[0] = r_x;
    this->r[1] = r_y;
    this->r[2] = r_z;
    
    
    double dr_x = (dO[0] * (cos(this->omega) * cos(this->big_omega)
                       - sin(this->omega) * sin(this->big_omega) 
                       * cos(this->i))
                - dO[1] * (sin(this->omega) * cos(this->big_omega)
                       + cos(this->omega) * cos(this->i) 
                       * sin(this->big_omega)));
    double dr_y = (dO[0] * (cos(this->omega) * sin(this->big_omega)
                       + sin(this->omega) * cos(this->big_omega) 
                       * cos(this->i))
                + dO[1] * (cos(this->omega) * cos(this->big_omega)
                       * cos(this->i) - sin(this->big_omega) 
                       * sin(this->omega)));
    double dr_z = (dO[0] * sin(this->omega) * sin(this->i)
                 + dO[1] * cos(this->omega) * sin(this->i));
    this->dr[0] = dr_x;
    this->dr[1] = dr_y;
    this->dr[2] = dr_z;
}

double* BodyPosition::getPosition() {

    return this->r;
}

double* BodyPosition::getVelocity() {
    
    return this->dr;
}

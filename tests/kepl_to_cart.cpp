#include "kepl_to_cart.h"
#include <iostream>
#include <cmath>
#include <tuple>

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

BodyPosition::BodyPosition(double* r, double* dr, double mu)
{
    this->r[0] = r[0];
    this->r[1] = r[1];
    this->r[2] = r[2];
    this->dr[0] = dr[0];
    this->dr[1] = dr[1];
    this->dr[2] = dr[2];
    
    this->mu = mu;
}

void BodyPosition::forwardConvertation() {
    double E  = this->m0 + this->e * sin(this->m0)/(1 - this->e * cos(this->m0));
    double x  = sqrt(1 - this->e) * cos(E/2);
    double y  = sqrt(1 + this->e) * sin(E/2);
    double nu = 2 * atan2(y, x); 

    double rc = this->a * (1 - this->e * cos(E));
    
    double O[3]     = {rc * cos(nu), rc * sin(nu), 0};
    double param    = (sqrt(this->mu * this->a))/rc;
    double dO_x     = - sin(E);
    double dO_y     = sqrt(1 - pow(this->e, 2)) * cos(E);
    double dO[3]    = {param * dO_x, param * dO_y, 0};

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
    double r_z = (O[0]  * sin(this->omega) * sin(this->i)
                + O[1]  * cos(this->omega) * sin(this->i));
    this->r[0] = r_x;
    this->r[1] = r_y;
    this->r[2] = r_z;
    
    
    double dr_x = (dO[0] * (cos(this->omega) * cos(this->big_omega)
                          - sin(this->omega) * sin(this->big_omega) 
                          * cos(this->i))
                - dO[1]  * (sin(this->omega) * cos(this->big_omega)
                          + cos(this->omega) * cos(this->i) 
                          * sin(this->big_omega)));
    double dr_y = (dO[0] * (cos(this->omega) * sin(this->big_omega)
                          + sin(this->omega) * cos(this->big_omega) 
                          * cos(this->i))
                + dO[1]  * (cos(this->omega) * cos(this->big_omega)
                          * cos(this->i) - sin(this->big_omega) 
                          * sin(this->omega)));
    double dr_z = (dO[0]  * sin(this->omega) * sin(this->i)
                 + dO[1]  * cos(this->omega) * sin(this->i));
    this->dr[0] = dr_x;
    this->dr[1] = dr_y;
    this->dr[2] = dr_z;
}

std::tuple<double, double, double> BodyPosition::vectorProduct(double* a, double* b) {

    double c[3] = {a[1] * b[2] - a[2] * b[1],
                  -a[0] * b[2] + a[2] * b[0],
                   a[0] * b[1] - a[1] * b[0]};
    return std::make_tuple(c[0], c[1], c[2]);
}

double BodyPosition::vectorMagnitude(double* a) {

    double mod_a = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    return mod_a;
}

double BodyPosition::scalarProduct(double* a, double* b) {

    double c = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    return c;
}

void BodyPosition::backConvertation() {

    auto h0 = vectorProduct(this->r, this->dr);
    double h[3] = {std::get<0>(h0), std::get<1>(h0), std::get<2>(h0)};
    //auto e0 = vectorProduct(this->r, h);
    double mod_r = vectorMagnitude(this->r);
    double mod_dr = vectorMagnitude(this->dr);
    //double e[3] = {std::get<0>(e0)/this->mu - this->r[0]/mod_r,
    //               std::get<1>(e0)/this->mu - this->r[1]/mod_r,
    //               std::get<2>(e0)/this->mu - this->r[2]/mod_r};
    double par_1 = mod_dr*mod_dr/this->mu - 1/mod_r;
    double par_2 = - scalarProduct(this->r, this->dr)/this->mu;
    double e[3]  = {par_1 * this->r[0] + par_2 * this->dr[0],
                    par_1 * this->r[1] + par_2 * this->dr[1],
                    par_1 * this->r[2] + par_2 * this->dr[2]};
                    
    this->e = vectorMagnitude(e);
    std::cout << this->e << std::endl;
    this->i = acos(h[2]/vectorMagnitude(h));
    
    double nu = 0;
    
    if (scalarProduct(this->r, this->dr) >= 0)
    {   
        nu = acos(scalarProduct(e, this->r)/(this->e * mod_r));
    }
    else
    {
        nu = 2*M_PI - acos(scalarProduct(e, this->r)/(this->e * mod_r));
    }

    double param = sqrt((1 + this->e)/(1 - this->e));
    //double E = 2 * atan(tan(nu/2)/param);
    this->a = 1/(2/mod_r - (mod_dr*mod_dr)/this->mu);
    double E = acos((this->a - mod_r)/(this->a*this->e));
    //double E = 2 * atan(tan(nu/2)/(sqrt((1+this->e)/(1-this->e))));
    this->m0 = E - (this->e * sin(E));

    double n[3] = {-h[1], h[0], 0};
    double mod_n = vectorMagnitude(n);

    if (n[1] >= 0)
    {
        this->big_omega = acos(n[0]/mod_n);
    }
    else
    {
        this->big_omega = 2*M_PI - acos(n[0]/mod_n);
    }

    if (e[2] >= 0)
    {
        this->omega = acos(scalarProduct(n, e)/(mod_n * this->e));
    }
    else
    {
        this->omega = 2*M_PI - acos(scalarProduct(n, e)/(mod_n * this->e));
    }
}

double* BodyPosition::getPosition() {

    return this->r;
}

double* BodyPosition::getVelocity() {
    
    return this->dr;
}

double BodyPosition::getA() {

    return this->a;
}

double BodyPosition::getE() {

    return this->e;
}

double BodyPosition::getOmega() {

    return this->omega;
}

double BodyPosition::getBigOmega() {

    return this->big_omega;
}

double BodyPosition::getI() {

    return this->i;
}

double BodyPosition::getM0() {

    return this->m0;
}

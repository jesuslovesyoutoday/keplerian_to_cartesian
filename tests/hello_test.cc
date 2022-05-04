#include <gtest/gtest.h>
#include "kepl_to_cart.h"
#include <iostream>
#include <cmath>

int main() {
    
    double a = 57909226541.52;
    double e = 0.20563593;
    double omega = 77.45779628 * M_PI/180; //1.35
    double big_omega = 48.33076593 * M_PI/180; // 0.84
    double i = 7.00497902 * M_PI/180; // 0.1221
    double m0 = 252.25032350 * M_PI/180; // 4.76
    //double mu = 3.986 * 10e14 ;
    double mu = 1.327124 * 10e20;
    //double mu = 1.98e+30;
    
    BodyPosition obj0 = BodyPosition(a, e, omega, big_omega, i, m0, mu);
    obj0.forwardConvertation();
    double* r  = obj0.getPosition();
    double* dr = obj0.getVelocity();
    
    double  r0[3] = {5.60950066e+10, 1.84707796e+10, -3.63983421e+09};
    double dr0[3] = {-2.30603729e+04, 4.04106010e+04, 5.41772383e+03};
    
    BodyPosition obj = BodyPosition(r, dr, mu);
    obj.backConvertation();
    double a_ = obj.getA();
    double e_ = obj.getE();
    double omega_ = obj.getOmega();
    double big_omega_ = obj.getBigOmega();
    double i_ = obj.getI();
    double m0_ = obj.getM0();
    
    EXPECT_EQ(round(a_*100)/100, round(a*100)/100);
    EXPECT_EQ(round(e_*100)/100, round(e*100)/100);
    EXPECT_EQ(round(omega_*100)/100, round(omega*100)/100);
    EXPECT_EQ(round(big_omega_*100)/100, round(big_omega*100)/100);
    EXPECT_EQ(round(i_*100)/100, round(i*100)/100);
    if (abs((360 - m0*180/M_PI)*(M_PI/180) - m0_) < 1)
    {
        EXPECT_EQ(round(m0_*100)/100, round((360 - m0*180/M_PI)*(M_PI/180)*100)/100);
    }
    else
    {
        EXPECT_EQ(round(m0_*100)/100, round(m0*100)/100);
    }
    return 0;
}

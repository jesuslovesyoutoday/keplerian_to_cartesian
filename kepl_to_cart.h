#ifndef KEPLERIAN_TO_CARTESIAN
#define KEPLERIAN_TO_CARTESIAN

#include <iostream>

class BodyPosition
{
    private:
            double a;          //< semi-major axis
            double e;          //< eccentricity
            double omega;      //< argument of periapsis
            double big_omega;  //< longitude of ascending node
            double i;          //< inclination
            double m0;         //< mean anomaly

            double r[3]  = {0, 0, 0};      //< position vector
            double dr[3] = {0, 0, 0};      //< velocity vector

            double mu;         //< gravitational parameter
    
    public:
            /*!
             * Constructor with Keplerian elements
             */
            BodyPosition(double a, double e, double omega, 
                         double big_omega, double i, double m0, double mu);
            /*!
             * Main translation function
             */
            void calculate();
            /*!
             * @return position vector
             */
            double* getPosition();
            /*!
             * @return velocity vector
             */
            double* getVelocity();
};

#endif // KEPLERIAN_TO_CARTESIAN

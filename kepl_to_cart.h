#ifndef KEPLERIAN_TO_CARTESIAN
#define KEPLERIAN_TO_CARTESIAN

#include <iostream>

class BodyPosition
{
    private:
            double a;
            double e;
            double omega;
            double big_omega;
            double i;
            double m0;         

            double r[3]  = {0, 0, 0};
            double dr[3] = {0, 0, 0};

            double mu;
    
    public:
            /*!
             * Constructor with Keplerian elements
             * (for forward convertation)
             * @param a         - semi-major axis, m
             * @param e         - eccentricity, 1
             * @param omega     - argument of periapsis, rad
             * @param big_omega - longitude of ascending node, rad
             * @param i         - inclination, rad
             * @param m0        - mean anomaly, rad
             * @param mu        - gravitational parameter, m^3*s^−2
             */
            BodyPosition(double a, double e, double omega, 
                         double big_omega, double i, double m0, double mu);
                         
            /*!
             * Constructor with cartesian coordinates
             * (for back convertation)
             * @param r  - position vector
             * @param dr - velocity vector
             * @param mu - gravitational parameter, m^3*s^−2
             */
            BodyPosition(double* r, double* dr, double mu);
            
            /*!
             * Func to translate from Keplerian elements
             * to cartesian coordinates
             */
            void forwardConvertation();
            
            /*!
             * Func for calculating vector product
             * @param a - first vector
             * @param b - second vector
             * @return    vector product
             */
            std::tuple<double, double, double> vectorProduct(double* a, double* b);
            
            /*!
             * Func for calculating magnitude
             * of given vector
             * @param a - given vector
             * @return    it's magnitude
             */
            double vectorMagnitude(double* a);
            
            /*!
             * Func for calculating scalar product
             * @param a - first vector
             * @param b - second vector
             * @return    scalar product
             */
            double scalarProduct(double* a, double* b);
            
            /*!
             * Func to translate from cartesian coords
             * to Keplerian elements
             */
            void backConvertation();
            
            /*!
             * @return position vector
             */
            double* getPosition();
            
            /*!
             * @return velocity vector
             */
            double* getVelocity();
            
            /*!
             * @return semi-major axis
             */
            double getA();
            
            /*!
             * @return eccentricity
             */
            double getE();
            
            /*!
             * @return argument of periapsis
             */
            double getOmega();
            
            /*!
             * @return longitude of ascending node
             */
            double getBigOmega();
            
            /*!
             * @return inclination
             */
            double getI();
            
            /*!
             * @return mean anomaly
             */
            double getM0();
};

#endif // KEPLERIAN_TO_CARTESIAN

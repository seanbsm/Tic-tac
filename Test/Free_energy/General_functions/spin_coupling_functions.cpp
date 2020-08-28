
#include "spin_coupling_functions.h"

double factorial(double n){
    if( n < 0 ) return -100.;

    return (n == 1. || n == 0.) ? 1. : factorial(n-1) * n ;
}

double CGcoeff(double J, double m, double J1, double m1, double J2, double m2){
    // (J1,m1) + (J2, m2) = (J, m)

    if( m != m1 + m2 ) return 0;

    double Jmin = abs(J1 - J2);
    double Jmax = J1+J2;

    if( J < Jmin || Jmax < J ) return 0;

    double s0 = (2*J+1.) * factorial(J+J1-J2) * factorial(J-J1+J2) * factorial(J1+J2-J) / factorial(J+J1+J2 + 1.);
    s0 = sqrt(s0);

    double s = factorial(J +m ) * factorial(J -m);
    double s1 = factorial(J1+m1) * factorial(J1-m1);
    double s2 = factorial(J2+m2) * factorial(J2-m2);
    s = sqrt(s * s1 * s2);

    //printf(" s0, s = %f , %f \n", s0, s);

    int kMax = std::min( std::min( J1+J2-J, J1 - m1), J2 + m2);

    double CG = 0.;
    for( int k = 0; k <= kMax; k++){
        double k1 = factorial(J1+J2-J-k);
        double k2 = factorial(J1-m1-k);
        double k3 = factorial(J2+m2-k);
        double k4 = factorial(J - J2 + m1 +k);
        double k5 = factorial(J - J1 - m2 +k);
        double temp = pow(-1, k) / (factorial(k) * k1 * k2 * k3 * k4 * k5);
        if( k1 == -100. || k2 == -100. || k3 == -100. || k4 == -100. || k5 == -100. ) temp = 0.;

        //printf(" %d | %f \n", k, temp);
        CG += temp;
    }

    return s0*s*CG;
}

double ThreeJSymbol(double J1, double m1, double J2, double m2, double J3, double m3){
    // ( J1 J2 J3 ) = (-1)^(J1-J2 - m3)/ sqrt(2*J3+1) * CGcoeff(J3, -m3, J1, m1, J2, m2) 
    // ( m1 m2 m3 )

    return pow(-1, J1 - J2 - m3)/sqrt(2*J3+1) * CGcoeff(J3, -m3, J1, m1, J2, m2);
}

double SixJSymbol(double J1, double J2, double J3, double J4, double J5, double J6){
    // coupling of j1, j2, j3 to J-J1
    // J1 = j1
    // J2 = j2
    // J3 = j12 = j1 + j2
    // J4 = j3
    // J5 = J = j1 + j2 + j3
    // J6 = j23 = j2 + j3

    //check triangle condition
    if( J3 < abs(J1 - J2 ) || J1 + J2 < J3 ) return 0; 
    if( J6 < abs(J2 - J4 ) || J2 + J4 < J6 ) return 0;
    if( J5 < abs(J1 - J6 ) || J1 + J6 < J5 ) return 0;
    if( J5 < abs(J3 - J4 ) || J3 + J4 < J5 ) return 0;

    double sixJ = 0;

    for( float m1 = -J1; m1 <= J1 ; m1 = m1 + 1){
        for( float m2 = -J2; m2 <= J2 ; m2 = m2 + 1){
            for( float m3 = -J3; m3 <= J3 ; m3 = m3 + 1){
                for( float m4 = -J4; m4 <= J4 ; m4 = m4 + 1){
                    for( float m5 = -J5; m5 <= J5 ; m5 = m5 + 1){
                        for( float m6 = -J6; m6 <= J6 ; m6 = m6 + 1){
                            double f = (J1 - m1) + (J2 - m2) + (J3 - m3) + (J4 - m4) + (J5 - m5) + (J6 - m6);

                            double a1 = ThreeJSymbol( J1, -m1, J2, -m2, J3, -m3); // J3 = j12 
                            double a2 = ThreeJSymbol( J1, m1, J5, -m5, J6, m6); // J5 = j1 + (J6 = j23)
                            double a3 = ThreeJSymbol( J4, m4, J2, m2, J6, -m6); // J6 = j23
                            double a4 = ThreeJSymbol( J4, -m4, J5, m5, J3, m3); // J5 = j3 + j12

                            double a = a1 * a2 * a3 * a4;
                            //if( a != 0 ) printf("%4.1f %4.1f %4.1f %4.1f %4.1f %4.1f | %f \n", m1, m2, m3, m4, m5, m6, a);

                            sixJ += pow(-1, f) * a1 * a2 * a3 * a4;
                        }
                    }
                }
            }
        }
    }

    return sixJ;
}


double NineJSymbol( double J1, double J2, double J3, double J4, double J5, double J6, double J7, double J8, double J9){
    double gMin = std::min( std::min (std::min( abs(J1 - J2 ), abs(J4 - J5)) , abs( J4 - J6 )) , abs(J7 - J8));
    double gMax = std::max( std::max (std::max( J1+J2, J4+J5), J3+J6), J7+J8);

    //printf(" gMin, gMax = %f %f \n", gMin, gMax);

    double nineJ = 0;
    for (float g = gMin; g <= gMax ; g = g + 1){
        double f = pow(-1, 2*g) * (2*g+1);

        double s1 = SixJSymbol(J1, J4, J7, J8, J9, g);
        if (s1 == 0) continue;

        double s2 = SixJSymbol(J2, J5, J8, J4, g, J6);
        if (s2 == 0) continue;

        double s3 = SixJSymbol(J3, J6, J9, g, J1, J2);
        if (s3 == 0) continue;

        nineJ += f* s1*s2*s3;
    }

    return nineJ;
}
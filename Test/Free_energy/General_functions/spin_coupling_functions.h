#ifndef SPIN_COUPLING_FUNCTIONS_H
#define SPIN_COUPLING_FUNCTIONS_H

#include <stdlib.h>
#include <cmath>
#include <algorithm>

double factorial(double n);
double CGcoeff(double J, double m, double J1, double m1, double J2, double m2);
double ThreeJSymbol(double J1, double m1, double J2, double m2, double J3, double m3);
double SixJSymbol(double J1, double J2, double J3, double J4, double J5, double J6);
double NineJSymbol( double J1, double J2, double J3, double J4, double J5, double J6, double J7, double J8, double J9);

#endif // SPIN_COUPLING_FUNCTIONS_H

#include "coupling_coefficients.h"

/* Clebsch-Gordan coefficients (it simply calls the Wigner-3j function and converts) */
double clebsch_gordan(int two_j1, int two_j2, int two_j3,
					  int two_m1, int two_m2, int two_m3){
	double result = 0;
	if ((abs(two_m1) <= two_j1) && (abs(two_m2) <= two_j2) && (abs(two_m1) <= two_j1))
	{
		result = gsl_sf_pow_int(-1, (two_j1 - two_j2 + two_m3) / 2) * sqrt(two_j3 + 1) * gsl_sf_coupling_3j (two_j1, two_j2, two_j3, two_m1, two_m2, -two_m3);
	}
	return result;
}

/* Wigner-9j symbols */
double wigner_3j(int two_j1, int two_j2, int two_j3,
				 int two_m1, int two_m2, int two_m3){
	return gsl_sf_coupling_3j (two_j1, two_j2, two_j3, two_m1, two_m2, two_m3);
}

/* Wigner-9j symbols */
double wigner_6j(int two_j1, int two_j2, int two_j3,
				 int two_j4, int two_j5, int two_j6){
	return gsl_sf_coupling_6j(two_j1, two_j2, two_j3, two_j4, two_j5, two_j6);
}

/* Wigner-9j symbols */
double wigner_9j(int two_j1, int two_j2, int two_j3,
				 int two_j4, int two_j5, int two_j6,
				 int two_j7, int two_j8, int two_j9){
	return gsl_sf_coupling_9j(two_j1, two_j2, two_j3, two_j4, two_j5, two_j6, two_j7, two_j8, two_j9);
}
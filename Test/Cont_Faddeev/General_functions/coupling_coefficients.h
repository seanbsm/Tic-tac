
#include <math.h>
#include <stdlib.h>

/* GSL functionality - contains efficient functions for coupling coefficients */
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf.h>

double clebsch_gordan(int two_j1, int two_j2, int two_j3,
                      int two_m1, int two_m2, int two_m3);

double wigner_3j(int two_j1, int two_j2, int two_j3,
                 int two_m1, int two_m2, int two_m3);

double wigner_6j(int two_j1, int two_j2, int two_j3,
                 int two_j4, int two_j5, int two_j6);

double wigner_9j(int two_j1, int two_j2, int two_j3,
                 int two_j4, int two_j5, int two_j6,
                 int two_j7, int two_j8, int two_j9);
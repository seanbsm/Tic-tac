#ifndef FADDEEV_ITERATOR_H
#define FADDEEV_ITERATOR_H

#include <math.h>

#include "General_functions/coupling_coefficients.h"

double p1(double q, double qp, double x);
double p2(double q, double qp, double x);

void iterate_faddeev(double** state_3N_symm_array,
                     int& Np, double** p_array, double** wp_array,
                     int& Nq, double** q_array, double** wq_array,
                     int& Nalpha, 
                     int** L_2N_array,
                     int** S_2N_array,
                     int** J_2N_array,
                     int** T_2N_array,
                     int** l_3N_array,
                     int** two_j_3N_array);

#endif // FADDEEV_ITERATOR_H
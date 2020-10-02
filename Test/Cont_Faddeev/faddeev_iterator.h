#ifndef FADDEEV_ITERATOR_H
#define FADDEEV_ITERATOR_H

#include <math.h>

double p1(double q, double qp, double x);
double p2(double q, double qp, double x);

void iterate_faddeev(double** state_3N_symm_array,
                     int& Np, double** p, double** wp,
                     int& Nq, double** q, double** wq,
                     int& Nalpha, int** L_2N, int** S_2N, int** J_2N, int** T_2N, int** l_3N, int** two_j_3N);

#endif // FADDEEV_ITERATOR_H
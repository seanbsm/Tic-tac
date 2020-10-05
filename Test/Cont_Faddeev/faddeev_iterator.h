#ifndef FADDEEV_ITERATOR_H
#define FADDEEV_ITERATOR_H

#include <math.h>

#include "lippmann_schwinger_solver.h"
#include "General_functions/coupling_coefficients.h"

void modified_gram_schmidt(double* state_matrix, double* state_basis, int N);

int kronecker_delta(int i, int j);

void S_spline (double* S, double* p_par, int N_par);

void iterate_faddeev(double* state_3N_symm_array,
                     double* P123_array,
                     int Np, double* p_array, double* wp_array,
                     int Nq, double* q_array, double* wq_array,
                     int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array, int* l_3N_array, int* two_j_3N_array,
                     int two_T, int two_J, int PAR);

#endif // FADDEEV_ITERATOR_H
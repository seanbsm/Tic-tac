#ifndef FADDEEV_ITERATOR_H
#define FADDEEV_ITERATOR_H

#include <math.h>

#include "lippmann_schwinger_solver.h"
#include "auxiliary.h"
#include "General_functions/coupling_coefficients.h"
#include "Interactions/potential_model.h"

void modified_gram_schmidt(double* state_matrix, double* state_basis, int N);

void iterate_faddeev(double* state_3N_symm_array,
                     int Np, double* p_array, double* wp_array,
                     int Nq, double* q_array, double* wq_array,
                     int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array, int* l_3N_array, int* two_j_3N_array,
                     int idx_alpha_proj, int idx_p_proj, int idx_q_proj,
                     int two_T, int two_J, int PAR,
                     potential_model* pot_ptr_np,
                     potential_model* pot_ptr_nn);

#endif // FADDEEV_ITERATOR_H
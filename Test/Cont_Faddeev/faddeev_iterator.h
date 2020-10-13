#ifndef FADDEEV_ITERATOR_H
#define FADDEEV_ITERATOR_H

#include <math.h>
#include <complex>
#include "mkl.h"

#include "lippmann_schwinger_solver.h"
#include "auxiliary.h"
#include "General_functions/coupling_coefficients.h"
#include "Interactions/potential_model.h"

void find_eigenvalues(double* A, double* wr, int N);

void modified_gram_schmidt(double* state_matrix,
                           double* state_basis,
                           int num_states,
                           int num_state_elements);

void brute_force_lanczos_for_faddeev(double* states_array,
                                     double* physical_state_array,
                                     int& physical_state_idx,
                                     int num_states,
                                     int num_state_elements);

void make_G_array(double* G_array,
                  int Np, double* p_array,
                  int Nq, double* q_array,
                  int Nx, double* x_array,
                  int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array, int* l_3N_array, int* two_j_3N_array,
                  int two_T, int two_J);

double interpolate_using_spline_array(double* S_array, int Np, double* p_array, int idx_pi, double p, int idx_p);

void calculate_angular_quadrature_grids(double* x_array, double* wx_array, int Nx);

void iterate_faddeev(double* state_3N_symm_array,
                     double* G_array,
                     int Np, double* p_array, double* wp_array,
                     int Nq, double* q_array, double* wq_array,
                     int Nx, double* x_array, double* wx_array,
                     int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array, int* l_3N_array, int* two_j_3N_array,
                     int idx_alpha, int idx_p, int idx_q,
                     int two_T, int two_J, int PAR,
                     potential_model* pot_ptr_nn,
                     potential_model* pot_ptr_np);

#endif // FADDEEV_ITERATOR_H
#ifndef FADDEEV_ITERATOR_H
#define FADDEEV_ITERATOR_H

#include <math.h>
#include <complex>
#include "mkl.h"

#include "error_management.h"
#include "lippmann_schwinger_solver.h"
#include "auxiliary.h"
#include "General_functions/coupling_coefficients.h"
#include "Interactions/potential_model.h"

void find_eigenvalues(double* A, double* wr, double* vr, int N);

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

double extract_potential_element_from_array(int L, int Lp, int J, int S, bool coupled, double* V_array);

void calculate_potential_matrices_array(double* V_mat_array,
                                        int Np, double* p_array, double* wp_array,
                                        int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array,
                                        potential_model* pot_ptr_nn,
                                        potential_model* pot_ptr_np);

void calculate_faddeev_convergence(double* state_array,
                                   double* P123_array,
                                   int Np, double* p_array, double* wp_array,
                                   int Nq, double* q_array, double* wq_array,
                                   int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array, int* l_3N_array, int* two_j_3N_array,
                                   int two_T, int two_J, int PAR, int J_2N_max,
                                   potential_model* pot_ptr_nn,
                                   potential_model* pot_ptr_np);

double iterate_faddeev(double Z,
                       double* P123_array,
                       double* V_unco_array,
                       double* V_coup_array,
                       int Np, double* p_array, double* wp_array,
                       int Nq, double* q_array, double* wq_array,
                       int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array, int* l_3N_array, int* two_j_3N_array,
                       int idx_alpha_r, int idx_p_r, int idx_q_r,
                       int idx_alpha_c, int idx_p_c, int idx_q_c,
                       int two_T, int two_J, int PAR,
                       potential_model* pot_ptr_nn,
                       potential_model* pot_ptr_np);

#endif // FADDEEV_ITERATOR_H
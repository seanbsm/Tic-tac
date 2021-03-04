#ifndef MAKE_SWP_STATES_H
#define MAKE_SWP_STATES_H

#include <iostream>
#include <math.h>
#include <complex>
#include "mkl.h"

#include "constants.h"
#include "error_management.h"

void diagonalize_real_symm_matrix(float *A, float *w, float *z, int N);
void diagonalize_real_symm_matrix(double *A, double *w, double *z, int N);

void construct_free_hamiltonian(double* H0_WP_array,
								int Np_WP, double* p_WP_array);

/* Constructs a NN-pair full Hamiltonians as upper-triangular from given arrays */
void construct_full_hamiltonian(double* mat_ptr_H,
								double* mat_ptr_H0,
								double* mat_ptr_V,
								int     mat_dim);

/* This function reorders the eigenvalues and corresponding vectors to correspond to coupled channels.
 * To understand why it's best to read up on WPCD theory - it's too long to explain here. */
void reorder_coupled_eigenspectrum(double* eigenvalues,
								   double* eigenvectors,
								   int     Np_WP);

void look_for_unphysical_bound_states(double* eigenvalues,
						   			  int     mat_dim,
						   			  bool    chn_3S1);

void make_swp_bin_boundaries(double* eigenvalues,
							 double* e_SWP_array,
							 int	 Np_WP,
							 bool    coupled,
							 bool    chn_3S1);

void make_swp_states(double* e_SWP_unco_array,
					 double* e_SWP_coup_array,
					 double* C_WP_unco_array,
					 double* C_WP_coup_array,
					 double* V_WP_unco_array,
					 double* V_WP_coup_array,
					 int Np_WP, double* p_WP_array,
					 int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array,
					 int J_2N_max);

#endif // MAKE_SWP_STATES_H
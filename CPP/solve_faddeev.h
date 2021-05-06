#ifndef SOLVE_FADDEEV_H
#define SOLVE_FADDEEV_H

#include <iostream>
#include <stdlib.h>     /* div, div_t */
#include <math.h>
#include <complex>
#include <omp.h>

/* Time-keeping modules */
#include <chrono>
#include <ctime>

#include "constants.h"
#include "type_defs.h"
#include "store_functionality.h"
#include "error_management.h"
#include "General_functions/matrix_routines.h"

void solve_faddeev_equations(cdouble*  U_array,
							 cdouble*  G_array,
							 double*   P123_sparse_val_array,
							 int*      P123_sparse_row_array,
							 int*      P123_sparse_col_array,
							 size_t    P123_sparse_dim,
							 double*   C_WP_unco_array,
							 double*   C_WP_coup_array,
							 double*   V_WP_unco_array,
							 double*   V_WP_coup_array,
							 int*	   q_com_idx_array,	   size_t num_q_com,
					  		 int*      deuteron_idx_array, size_t num_deuteron_states,
							 int       J_2N_max,
							 size_t    Nq_WP,
							 size_t    Np_WP,
							 size_t    Nalpha,
							 int*      L_2N_array,
							 int*      S_2N_array,
							 int*      J_2N_array,
							 int*      T_2N_array,
							 int*      L_1N_array, 
							 int*      two_J_1N_array);

/* Create array of pointers to C^T matrices for product (C^T)PVC in row-major format */
void create_CT_row_maj_3N_pointer_array(double** CT_RM_array,
										double*  C_WP_unco_array,
										double*  C_WP_coup_array,
										size_t   Np_WP,
										int      J_2N_max,
										size_t   Nalpha,
										int*     L_2N_array,
										int*     S_2N_array,
										int*     J_2N_array,
										int*     T_2N_array,
									    int*	 L_1N_array,
									    int*	 two_J_1N_array);

/* Restructures NN coupled matrix as 4 seperate matrices
 * !!! WARNING: COLUMN-MAJOR ALGORITHM !!! */
void restructure_coupled_VC_product(double* VC_product, size_t Np_WP);

/* Create array of pointers to VC-product matrices for product (C^T)PVC in column-major format*/
void create_VC_col_maj_3N_pointer_array(double** VC_CM_array,
										double*  C_WP_unco_array,
										double*  C_WP_coup_array,
										double*  V_WP_unco_array,
										double*  V_WP_coup_array,
										size_t   Np_WP,
										int      J_2N_max,
										size_t   Nalpha,
										int*     L_2N_array,
										int*     S_2N_array,
										int*     J_2N_array,
										int*     T_2N_array,
									    int*	 L_1N_array,
									    int*	 two_J_1N_array);

void PVC_col_brute_force(double*  col_array,
					     size_t   idx_alpha_c, size_t idx_p_c, size_t idx_q_c,
					     size_t   Nalpha,      size_t Nq_WP,   size_t Np_WP,
					     double** VC_CM_array,
					     double*  P123_val_array,
					     size_t*  P123_row_array,
					     int*     P123_col_array,
					     size_t   P123_dim);

void CPVC_col_brute_force(double*  col_array,
						  size_t   idx_alpha_c, size_t idx_p_c, size_t idx_q_c,
						  size_t   Nalpha,      size_t Nq_WP,   size_t Np_WP,
						  double** CT_RM_array,
						  double** VC_CM_array,
						  double*  P123_val_array,
						  size_t*  P123_row_array,
						  int*     P123_col_array,
						  size_t   P123_dim);

void PVC_col_calc_test(size_t   Nalpha,
					   size_t 	Nq_WP,
					   size_t 	Np_WP,
					   double** VC_CM_array,
					   double*  P123_sparse_val_array,
					   int*     P123_sparse_row_array,
					   size_t*  P123_sparse_col_array,
					   size_t   P123_sparse_dim);

void CPVC_col_calc_test(size_t   Nalpha,
						size_t 	 Nq_WP,
						size_t 	 Np_WP,
						double** CT_RM_array,
						double** VC_CM_array,
						double*  P123_sparse_val_array,
						int*     P123_sparse_row_array,
						size_t*  P123_sparse_col_array,
						size_t   P123_sparse_dim);

#endif // SOLVE_FADDEEV_H
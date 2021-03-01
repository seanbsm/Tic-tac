#ifndef SOLVE_FADDEEV_H
#define SOLVE_FADDEEV_H

#include <iostream>
#include <stdlib.h>     /* div, div_t */
#include <math.h>
#include <complex>

/* Time-keeping modules */
#include <chrono>
#include <ctime>

#include "constants.h"
#include "type_defs.h"
#include "error_management.h"
#include "General_functions/matrix_routines.h"

void pade_method_solve();
void direct_sparse_solve();

/* Solves the Faddeev equations
 * U = P*V + P*V*G*U
 * on the form L*U = R, where L and R are the left-
 * and right-handed sides of the equations, given by 
 * L = 1 - P*V*G
 * R = P*V
 * Since G is expressed in an SWP basis, we also must include the basis-transormation matrices C */
void direct_dense_solve(cdouble* U_array,
						cdouble* G_array,
						double* P123_array,
						double* C_WP_unco_array,
						double* C_WP_coup_array,
						double* V_WP_unco_array,
						double* V_WP_coup_array,
						int Nq_WP,
						int Np_WP,
						int Nalpha,
						int* L_2N_array,
						int* S_2N_array,
						int* J_2N_array,
						int* T_2N_array,
						int* L_1N_array, 
						int* two_J_1N_array,
						int* two_J_3N_array,
						int* two_T_3N_array);

void solve_faddeev_equations(cdouble*  U_array,
							 cdouble*  G_array,
							 double*   P123_sparse_val_array,
							 int*      P123_sparse_row_array,
							 int*      P123_sparse_col_array,
							 int       P123_sparse_dim,
							 double*   C_WP_unco_array,
							 double*   C_WP_coup_array,
							 double*   V_WP_unco_array,
							 double*   V_WP_coup_array,
							 int       J_2N_max,
							 int       Nq_WP,
							 int       Np_WP,
							 int       Nalpha,
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
										int      Np_WP,
										int      J_2N_max,
										int      Nalpha,
										int*     L_2N_array,
										int*     S_2N_array,
										int*     J_2N_array,
										int*     T_2N_array);

/* Restructures NN coupled matrix as 4 seperate matrices
 * !!! WARNING: COLUMN-MAJOR ALGORITHM !!! */
void restructure_coupled_VC_product(double* VC_product, int Np_WP);

/* Create array of pointers to VC-product matrices for product (C^T)PVC in column-major format*/
void create_VC_col_maj_3N_pointer_array(double** VC_CM_array,
										double*  C_WP_unco_array,
										double*  C_WP_coup_array,
										double*  V_WP_unco_array,
										double*  V_WP_coup_array,
										int      Np_WP,
										int      J_2N_max,
										int      Nalpha,
										int*     L_2N_array,
										int*     S_2N_array,
										int*     J_2N_array,
										int*     T_2N_array);

#endif // SOLVE_FADDEEV_H
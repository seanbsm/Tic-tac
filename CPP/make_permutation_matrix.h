#ifndef PERMUTATION_OPERATORS_H
#define PERMUTATION_OPERATORS_H

#include <string>
#include <cstring> //std::strcpy
#include <iostream>
#include <fstream>
#include <complex>
#include <numeric>
#include <algorithm>
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream
#include <vector>

/* Modules used in benchmark code */
#include "mkl.h"
#include <math.h>
#include "auxiliary.h"

/* GSL functionality */
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf.h>

/* Time-keeping modules */
#include <chrono>
#include <ctime>

/* Project headers */
#include "error_management.h"
#include "store_functionality.h"
#include "General_functions/matrix_routines.h"

void calculate_permutation_matrix_for_3N_channel(double** P123_val_dense_array,
												 double** P123_val_sparse_array,
												 int**    P123_row_array,
												 int**    P123_col_array,
												 size_t&  P123_dim,
												 bool     use_dense_format,
												 bool     production_run,
												 int      Np_WP, double *p_array_WP_bounds,
												 int      Nq_WP, double *q_array_WP_bounds,
												 int      Nx, double* x_array, double* wx_array,
												 int      Nphi,
												 int      J_2N_max,
												 int      Nalpha,
												 int*     L_2N_array,
												 int*     S_2N_array,
												 int*     J_2N_array,
												 int*     T_2N_array,
												 int*     L_1N_array,
												 int*     two_J_1N_array,
												 int      two_J_3N,
												 int      two_T_3N,
												 int      P_3N);

void calculate_permutation_matrices_for_all_3N_channels(double** P123_sparse_val_array,
														int**    P123_sparse_row_array,
														int**    P123_sparse_col_array,
														size_t&  P123_sparse_dim_array,
														bool     production_run,
														int      Np_WP, double *p_array_WP_bounds,
														int      Nq_WP, double *q_array_WP_bounds,
														int      Nx, double* x_array, double* wx_array,
														int      Nphi,
														int      J_2N_max,
														int      Nalpha,
														int*     L_2N_array,
														int*     S_2N_array,
														int*     J_2N_array,
														int*     T_2N_array,
														int*     L_1N_array,
														int*     two_J_1N_array,
														int      two_J_3N,
														int      two_T_3N,
												 		int      P_3N);

#endif // PERMUTATION_OPERATORS_H
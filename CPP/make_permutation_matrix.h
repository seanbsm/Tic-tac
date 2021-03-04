#ifndef PERMUTATION_OPERATORS_H
#define PERMUTATION_OPERATORS_H

#include <string>
#include <cstring> //std::strcpy
#include <iostream>
#include <fstream>
#include <complex>
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

/* .h5-file functionality */
#include "hdf5/serial/hdf5.h"
#include "hdf5/serial/hdf5_hl.h"

/* Project headers */
#include "error_management.h"
#include "General_functions/matrix_routines.h"

typedef struct pw_table
{
	int alpha_pwtable;
	int L_pwtable;
	int S_pwtable;
	int J_pwtable;
	int T_pwtable;
	int l_pwtable;
	int twoj_pwtable;
	int twoJtotal_pwtable;
	int PARtotal_pwtable;
	int twoTtotal_pwtable;
} pw_table;

typedef struct alpha_table
{
	int alpha_alphaprime_index;
	int alpha_index;
	int alphaprime_index;
} alpha_table;


typedef struct pmesh_table
{
	int index_ptable;
	double mesh_ptable;
	double weight_ptable;
} pmesh_table;


typedef struct qmesh_table
{
	int index_qtable;
	double mesh_qtable;
	double weight_qtable;
} qmesh_table;

void check_h5_read_call(herr_t ret);
void check_h5_read_table_call(herr_t ret);
void check_h5_close_call(herr_t ret);

void read_P123_h5_data_file(std::string file_path,
							double* P123,
							int Np_3N, double* p_3N_3, double* wp_3N_3,
							int Nq_3N, double* q_3N_3, double* wq_3N_3);

void get_h5_P123_dimensions(std::string file_path, int& Nalpha, int& Np, int& Nq);

void calculate_permutation_matrix_for_3N_channel(double** P123_val_dense_array,
												 double** P123_val_sparse_array,
												 int** P123_row_array,
												 int** P123_col_array,
												 int& P123_dim,
												 bool use_dense_format,
												 int Nq, double* q_array, double* wq_array, int Np_per_WP, int Np_WP, double *p_array_WP_bounds,
												 int Np, double* p_array, double* wp_array, int Nq_per_WP, int Nq_WP, double *q_array_WP_bounds,
												 int Nx, double* x_array, double* wx_array,
												 int Nalpha,
												 int* L_2N_array,
												 int* S_2N_array,
												 int* J_2N_array,
												 int* T_2N_array,
												 int* L_1N_array,
												 int* two_J_1N_array,
												 int two_J_3N,
												 int two_T_3N);

void increase_sparse_array_size(double** sparse_array, int array_length, int sparse_step_length);
void increase_sparse_array_size(int**    sparse_array, int array_length, int sparse_step_length);

void reduce_sparse_array_size(double** sparse_array, int array_length, int sparse_dim);
void reduce_sparse_array_size(int**    sparse_array, int array_length, int sparse_dim);

void calculate_permutation_matrices_for_all_3N_channels(double** P123_sparse_val_array,
														int**    P123_sparse_row_array,
														int**    P123_sparse_col_array,
														int&     P123_sparse_dim_array,
														int  Nq, double* q_array, double* wq_array, int Np_per_WP, int Np_WP, double *p_array_WP_bounds,
														int  Np, double* p_array, double* wp_array, int Nq_per_WP, int Nq_WP, double *q_array_WP_bounds,
														int  Nx, double* x_array, double* wx_array,
														int  Nalpha,
														int* L_2N_array,
														int* S_2N_array,
														int* J_2N_array,
														int* T_2N_array,
														int* L_1N_array,
														int* two_J_1N_array,
														int  two_J_3N,
														int  two_T_3N);

#endif // PERMUTATION_OPERATORS_H
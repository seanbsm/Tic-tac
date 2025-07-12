#ifndef MATRIX_ROUTINES_H
#define MATRIX_ROUTINES_H

#include "../type_defs.h"
#include "templates.h"

#include <iostream>
#include <math.h>
#include <complex>
#include <vector>
#include <numeric>
#include <algorithm>

//#include "mkl.h"
//#include "mkl_dss.h"

//#include "cblas.h"


// extern "C" {
// 	#include <cblas.h>
// }

void simple_transpose_matrix_routine(double* mat_array, int mat_dim);

/* Converts a dense-format matrix to a sparse COO-format matrix */
void square_dense_to_sparse_COO_format_converter(int      mat_dim,
												 double*  mat_dense_array,
												 double** mat_sparse_val_array,
												 int**    mat_sparse_row_array,
												 int**    mat_sparse_col_array,
												 size_t&  mat_sparse_dim);

/* Converts a sparse COO-format matrix to a dense-format matrix */
void square_sparse_COO_to_dense_format_converter(int      mat_dim,
												 double** mat_dense_array,
												 double*  mat_sparse_val_array,
												 int*     mat_sparse_row_array,
												 int*     mat_sparse_col_array,
												 size_t   mat_sparse_dim);

/* Converts row-major sparse COO-format index arrays to CSR-format index arrays */
void coo_to_csr_format_converter(int*    idx_row_array_coo,
								 size_t* idx_row_array_csr,
								 size_t  mat_sparse_dim,
								 size_t  mat_dense_dim);
/* Converts column-major sparse COO-format index arrays to CSC-format index arrays */
void coo_to_csc_format_converter(int*    idx_col_array_coo,
								 size_t* idx_col_array_csc,
								 size_t  mat_sparse_dim,
								 size_t  mat_dense_dim);

/* Converts sparse random-order COO-format arrays to row-major COO-format */
void unsorted_sparse_to_coo_row_major_sorter(double** mat_val_array,
											 int**    mat_row_array,
											 int**    mat_col_array,
											 size_t   mat_sparse_dim,
											 int      mat_dense_dim);
/* Converts sparse random-order COO-format arrays to column-major COO-format */
void unsorted_sparse_to_coo_col_major_sorter(double** mat_val_array,
											 int**    mat_row_array,
											 int**    mat_col_array,
											 size_t   mat_sparse_dim,
											 int      mat_dense_dim);

/* Increase a full sparse array size of length array_length by sparse_step_length elements */
void increase_sparse_array_size(double** sparse_array, int array_length, int sparse_step_length);
void increase_sparse_array_size(int** sparse_array, int array_length, int sparse_step_length);

/* Decrease a non-full sparse array size of length array_length to sparse_dim elements */
void reduce_sparse_array_size(double** sparse_array, int array_length, size_t sparse_dim);
void reduce_sparse_array_size(int** sparse_array, int array_length, size_t sparse_dim);

#endif // MATRIX_ROUTINES_H
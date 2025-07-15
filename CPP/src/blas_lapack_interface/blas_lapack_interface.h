#ifndef BLAS_LAPACK_INTERFACE_H
#define BLAS_LAPACK_INTERFACE_H

// Local type definitions
typedef unsigned long long int tmp_int;

// Extra libraries
#include <complex>
#include <omp.h>

// BLAS complex float types (!!! CAN BE REMOVED?)
#define blas_complex_double std::complex<double>
#define blas_complex_float std::complex<float>
// LAPACK complex float types
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

/* -------------------------------------------------------------------------- */
/*                                BLAS functions                              */
/* -------------------------------------------------------------------------- */
void set_max_OMP_threads_for_BLAS();

void dot_MV(double *A, double *B, double *C, int N, int M);

std::complex<double> cdot_VV(blas_complex_float *X, blas_complex_float *Y, int N, int INCR_X, int INCR_Y);
std::complex<double> cdot_VV(blas_complex_double *X, blas_complex_double *Y, int N, int INCR_X, int INCR_Y);

/* Matrix-matrix dot-product with simplified input */
void dot_MM(float  *A, float  *B, float  *C, int N, int K, int M);
void dot_MM(double *A, double *B, double *C, int N, int K, int M);
/* Matrix-matrix dot-product with full input */
void dot_MM_full(bool       row_maj,
				 bool       no_trans_A,
				 bool       no_trans_B,
				 size_t		M,
				 size_t		N,
				 size_t		K,
				 float* 	A,
				 size_t 	lda,
				 float* 	B,
				 size_t 	ldb,
				 float* 	C,
				 size_t 	ldc);
void dot_MM_full(bool       row_maj,
				 bool       no_trans_A,
				 bool       no_trans_B,
				 size_t		M,
				 size_t		N,
				 size_t		K,
				 double* 	A,
				 size_t 	lda,
				 double* 	B,
				 size_t 	ldb,
				 double* 	C,
				 size_t 	ldc);
/* Matrix-matrix complex dot-product with simplified input */
void cdot_MM(blas_complex_float *A, blas_complex_float *B, blas_complex_float *C, int N, int K, int M);
void cdot_MM(blas_complex_double *A, blas_complex_double *B, blas_complex_double *C, int N, int K, int M);

/* Solves linear system A*X=B for dense NxN matrices, for matrix X */
void solve_MM(float *A, float *B, int N);
void solve_MM(double *A, double *B, int N);
void solve_MM(blas_complex_float *A, blas_complex_float *B, int N);
void solve_MM(blas_complex_double *A, blas_complex_double *B, int N);

// void dot_MM_sparse(double* A_v, long long int* A_c, long long int* A_r, double *B, double *C, int N, int K, int M, bool transpose_A);

void dimatcopy(
    bool    row_maj,
    bool    no_trans,
	int     rows,
	int     cols,
	double  alpha,
	double* A,
	int     lda,
	int     ldb);

/* -------------------------------------------------------------------------- */
/*                               LAPACK functions                             */
/* -------------------------------------------------------------------------- */

std::complex<double> determinant(std::complex<double>* A, int N);

/* Solves linear system A*X=B for sparse NxN matrix A and column B, for column X
 * The logic in this function follows exactly that listed in
 * https://scc.ustc.edu.cn/zlsc/sugon/intel/mkl/mkl_manual/GUID-24A41D50-4D0A-471B-9722-EA59E443ABC8.htm*/
void solve_MM_sparse(double*  A_val_array,
					 tmp_int* A_idx_row_array,
					 tmp_int* A_idx_col_array,
					 tmp_int  A_nnz,
					 double*  B_val_array,
					 tmp_int  mat_dim,
					 double*  sol_val_array);
void solve_MM_sparse(std::complex<double>*  A_val_array,
					 tmp_int* A_idx_row_array,
					 tmp_int* A_idx_col_array,
					 tmp_int  A_nnz,
					 std::complex<double>*  B_val_array,
					 tmp_int  mat_dim,
					 std::complex<double>*  sol_val_array);

#endif // DISK_IO_FUNCTIONALITY_H
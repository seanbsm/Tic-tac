
#include "blas_lapack_interface.h"

/* BLAS and LAPACK libraries */
#include <cblas.h>
#include <lapacke.h>


/* Use as many threads as possible in GEMM */
void set_max_OMP_threads_for_BLAS(){
	openblas_set_num_threads(omp_get_max_threads());
}

void dot_MV(double *A, double *B, double *C, int N, int M){
	double  beta  = 0.0;
	double  alpha = 1.0;
	//tmp_int incrx = 1.0;
	//tmp_int incry = 1.0;
	blasint incrx = 1.0;
	blasint incry = 1.0;
	char    trans = 'N';
	//tmp_int n = N;
	//tmp_int m = M;
	blasint n = N;
	blasint m = M;
	//cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, M, K, &alpha, A, K, B, M, &beta, C, M);
	//dgemv(&trans, &n, &m, &alpha, A, &m, B, &incrx, &beta, C, &incry);
	cblas_dgemv(CblasRowMajor, CblasNoTrans, n, m, alpha, A, m, B, incrx, beta, C, incry);
}


std::complex<double> cdot_VV(std::complex<float> *X, std::complex<float> *Y, int N, int INCR_X, int INCR_Y){
	std::complex<float> dot_product = 0;
	cblas_cdotu_sub(N, X, INCR_X, Y, INCR_Y, &dot_product);
	return dot_product;
}
std::complex<double> cdot_VV(std::complex<double> *X, std::complex<double> *Y, int N, int INCR_X, int INCR_Y){
	std::complex<double> dot_product = 0;
	cblas_zdotu_sub(N, X, INCR_X, Y, INCR_Y, &dot_product);
	return dot_product;
}

void dot_MM(float *A, float *B, float *C, int N, int K, int M){
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, M, K, 1.0, A, K, B, M, 0.0, C, M);
}
void dot_MM(double *A, double *B, double *C, int N, int K, int M){
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, M, K, 1.0, A, K, B, M, 0.0, C, M);
}

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
				 size_t 	ldc){
	/* Define layout */
	CBLAS_LAYOUT layout = CblasColMajor;
	if (row_maj){
		layout = CblasRowMajor;
	}
	/* Define transposition of matrix A */
	CBLAS_TRANSPOSE transA = CblasTrans;
	if (no_trans_A){
		transA = CblasNoTrans;
	}
	/* Define transposition of matrix B */
	CBLAS_TRANSPOSE transB = CblasTrans;	
	if (no_trans_B){
		transB = CblasNoTrans;
	}
	const float beta  = 0;
	const float alpha = 1;
	cblas_sgemm(layout, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
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
				 size_t 	ldc){
	/* Define layout */
	CBLAS_LAYOUT layout = CblasColMajor;
	if (row_maj){
		layout = CblasRowMajor;
	}
	/* Define transposition of matrix A */
	CBLAS_TRANSPOSE transA = CblasTrans;
	if (no_trans_A){
		transA = CblasNoTrans;
	}
	/* Define transposition of matrix B */
	CBLAS_TRANSPOSE transB = CblasTrans;	
	if (no_trans_B){
		transB = CblasNoTrans;
	}
	const double beta  = 0;
	const double alpha = 1;
	cblas_dgemm(layout, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

void cdot_MM(std::complex<float> *A, std::complex<float> *B, std::complex<float> *C, int N, int K, int M){
	std::complex<float> beta = {0,0};
	std::complex<float> alpha = {1,0};
	cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, M, K, &alpha, A, K, B, M, &beta, C, M);
}
void cdot_MM(std::complex<double> *A, std::complex<double> *B, std::complex<double> *C, int N, int K, int M){
	std::complex<double> beta = {0,0};
	std::complex<double> alpha = {1,0};
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, M, K, &alpha, A, K, B, M, &beta, C, M);
}


void dimatcopy(
    bool    row_maj,
    bool    no_trans,
	int     rows,
	int     cols,
	double  alpha,
	double* A,
	int     lda,
	int     ldb){
	/* Define layout */
	CBLAS_LAYOUT layout;
	if (row_maj){
		layout = CblasRowMajor;
	}
	else{
		layout = CblasColMajor;
	}
	/* Define transposition of matrix A */
	CBLAS_TRANSPOSE trans;
	if (no_trans){
		trans = CblasNoTrans;
	}
	else{
		trans = CblasTrans;
	}
	cblas_dimatcopy(layout, trans, rows, cols, alpha, A, lda, ldb);
}


/* Finds the eigenvalues and eigenvectors of a real, symmetric matrix A. For
 * simplicity here, A must be stored as an upper triangle. w will be filled with
 * the eigenvalues in ascending order, and z will be a matrix with the
 * corresponding eigenvectors (column by column). Lastly, we work with row-major
 * matrices, as is usual with C & C++ */
void diagonalize_real_symm_matrix(float *A, float *w, float *z, int N){
	char jobz = 'V';
	char uplo = 'U';
	
	LAPACKE_sspevd(LAPACK_ROW_MAJOR, jobz, uplo, N, A, w, z, N);
}
void diagonalize_real_symm_matrix(double *A, double *w, double *z, int N){
	char jobz = 'V';
	char uplo = 'U';
	
	LAPACKE_dspevd(LAPACK_ROW_MAJOR, jobz, uplo, N, A, w, z, N);
}
void solve_MM(float* A, float* B, int dim){
	char trans = 'N';
	lapack_int ipiv [dim];
	
	LAPACKE_sgetrf(LAPACK_ROW_MAJOR, dim, dim, A, dim, ipiv);
	LAPACKE_sgetrs(LAPACK_ROW_MAJOR, trans, dim, dim, A, dim, ipiv, B, dim);
}
void solve_MM(double* A, double* B, int dim){
	char trans = 'N';
	lapack_int ipiv [dim];
	
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, dim, dim, A, dim, ipiv);
	LAPACKE_dgetrs(LAPACK_ROW_MAJOR, trans, dim, dim, A, dim, ipiv, B, dim);
}
void solve_MM(std::complex<float> *A, std::complex<float> *B, int N){
	
	char trans = 'N';
	lapack_int ipiv [N];
	
	LAPACKE_cgetrf(LAPACK_ROW_MAJOR, N, N, A, N, ipiv);
	LAPACKE_cgetrs(LAPACK_ROW_MAJOR, trans, N, N, A, N, ipiv, B, N);
}
void solve_MM(std::complex<double> *A, std::complex<double> *B, int N){
	
	char trans = 'N';
	lapack_int ipiv [N];
	
	LAPACKE_zgetrf(LAPACK_ROW_MAJOR, N, N, A, N, ipiv);
	LAPACKE_zgetrs(LAPACK_ROW_MAJOR, trans, N, N, A, N, ipiv, B, N);
}

std::complex<double> determinant(std::complex<double>* A, int N){
	//char trans = 'N';
	lapack_int ipiv [N];
	
	/* Perform LU decomposition, A is overwritten by L */
	LAPACKE_zgetrf(LAPACK_ROW_MAJOR, N, N, A, N, ipiv);

	/* Diagonal of L equals determinant of A */
	std::complex<double> product_sum_diagonal = 1;
	for (int i=0; i<N; i++){        /* Figure out pivoting */
        if (ipiv[i]!=i+1){
            product_sum_diagonal *= -A[i*N+i];
        }
        else{
            product_sum_diagonal *= +A[i*N+i];
        }
	}

	return product_sum_diagonal;
}
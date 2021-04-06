
#include "matrix_routines.h"

void raise_error_DSS_MKL(_INTEGER_t error, const char* function_name){
    printf("Function %s returned error code %d\n", function_name, error);
	exit(1);
}

void dot_MM(float *A, float *B, float *C, int N, int K, int M){
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, M, K, 1.0, A, K, B, M, 0.0, C, M);
}
void dot_MM(double *A, double *B, double *C, int N, int K, int M){
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, M, K, 1.0, A, K, B, M, 0.0, C, M);
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

void solve_MM(float* A, float* B, int dim){
	char trans = 'N';
	long long int ipiv [dim];
	
	LAPACKE_sgetrf(LAPACK_ROW_MAJOR, dim, dim, A, dim, ipiv);
	LAPACKE_sgetrs(LAPACK_ROW_MAJOR, trans, dim, dim, A, dim, ipiv, B, dim);
}
void solve_MM(double* A, double* B, int dim){
	char trans = 'N';
	long long int ipiv [dim];
	
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, dim, dim, A, dim, ipiv);
	LAPACKE_dgetrs(LAPACK_ROW_MAJOR, trans, dim, dim, A, dim, ipiv, B, dim);
}
void solve_MM(std::complex<float> *A, std::complex<float> *B, int N){
	
	char trans = 'N';
	long long int ipiv [N];
	
	LAPACKE_cgetrf(LAPACK_ROW_MAJOR, N, N, A, N, ipiv);
	LAPACKE_cgetrs(LAPACK_ROW_MAJOR, trans, N, N, A, N, ipiv, B, N);
}
void solve_MM(std::complex<double> *A, std::complex<double> *B, int N){
	
	char trans = 'N';
	long long int ipiv [N];
	
	LAPACKE_zgetrf(LAPACK_ROW_MAJOR, N, N, A, N, ipiv);
	LAPACKE_zgetrs(LAPACK_ROW_MAJOR, trans, N, N, A, N, ipiv, B, N);
}

std::complex<double> determinant(std::complex<double>* A, int N){
	long long int ipiv [N];
	
	/* Perform LU decomposition, A is overwritten by L */
	LAPACKE_zgetrf(LAPACK_ROW_MAJOR, N, N, A, N, ipiv);

	/* Diagonal of L equals determinant of A */
	std::complex<double> product_sum_diagonal = 1;
	for (int i=0; i<N; i++){
		product_sum_diagonal *= A[i*N+i];
	}

	return product_sum_diagonal;
}

void solve_MM_sparse(double*  A_val_array,
					 MKL_INT* A_idx_row_array,
					 MKL_INT* A_idx_col_array,
					 MKL_INT  A_nnz,
					 double*  B_val_array,
					 MKL_INT  mat_dim,
					 double*  sol_val_array){
	
	bool print_content = false;

	_MKL_DSS_HANDLE_t dss_handle;
	_INTEGER_t        error;
	_CHARACTER_t statIn[] = "determinant";
	_DOUBLE_PRECISION_t statOut[5];

	MKL_INT  nRows     = mat_dim;
	MKL_INT  nCols     = mat_dim;
	MKL_INT  nNonZeros = A_nnz;
	MKL_INT  nRhs      = 1;
	MKL_INT* rowIndex  = A_idx_row_array;
	MKL_INT* columns   = A_idx_col_array;
	double*  values    = A_val_array;
	double*  rhs       = B_val_array;
	double*  solValues = sol_val_array;//[nRows];

	MKL_INT dss_opt  = MKL_DSS_MSG_LVL_WARNING + MKL_DSS_TERM_LVL_ERROR + MKL_DSS_ZERO_BASED_INDEXING;
	MKL_INT dss_sym  = MKL_DSS_NON_SYMMETRIC;
	MKL_INT dss_type = MKL_DSS_INDEFINITE;//MKL_DSS_POSITIVE_DEFINITE;

	/* Initialize DSS solver */
	error = dss_create(dss_handle,
					   dss_opt);
	if ( error != MKL_DSS_SUCCESS ){raise_error_DSS_MKL(error, "dss_handle");}
	if(print_content){std::cout << "Finished dss_create" << std::endl;}

	error = dss_define_structure(dss_handle,
								 dss_sym,
								 rowIndex,
								 nRows,
								 nCols,
								 columns,
								 nNonZeros);
	if ( error != MKL_DSS_SUCCESS ){raise_error_DSS_MKL(error, "dss_define_structure");}
	if(print_content){std::cout << "Finished dss_define_structure" << std::endl;}

	error = dss_reorder(dss_handle,
						dss_opt,
						0);
	if ( error != MKL_DSS_SUCCESS ){raise_error_DSS_MKL(error, "dss_reorder");}
	if(print_content){std::cout << "Finished dss_reorder" << std::endl;}

	error = dss_factor_real(dss_handle,
							dss_type,
							values );
	if ( error != MKL_DSS_SUCCESS ){raise_error_DSS_MKL(error, "dss_factor_real");}
	if(print_content){std::cout << "Finished dss_factor_real" << std::endl;}

	error = dss_solve_real(dss_handle,
						   dss_opt,
						   rhs,
						   nRhs,
						   solValues );
	if ( error != MKL_DSS_SUCCESS ){raise_error_DSS_MKL(error, "dss_solve_real");}
	if(print_content){std::cout << "Finished dss_solve_real" << std::endl;}

	if ( nRows < A_nnz ) {
		error = dss_statistics(dss_handle,
							   dss_opt,
							   statIn,
							   statOut);
		if ( error != MKL_DSS_SUCCESS ){raise_error_DSS_MKL(error, "dss_statistics");}
		
		if(print_content){
			printf(" - Determinant power is %g \n", statOut[0]);
			printf(" - Determinant base is %g \n", statOut[1]);
			printf(" - Determinant is %g \n", (pow(10.0,statOut[0]))*statOut[1]);
		}
	}
	if(print_content){std::cout << "Finished dss_statistics" << std::endl;}

	error = dss_delete(dss_handle,
					   dss_opt);
	if ( error != MKL_DSS_SUCCESS ){raise_error_DSS_MKL(error, "dss_delete");}
	if(print_content){std::cout << "Finished dss_delete" << std::endl;}

	if(print_content){
		printf(" Solution array: ");
		for(MKL_INT i = 0; i< nCols; i++)
			printf(" %g", solValues[i] );
		printf("\n");
	}
}
void solve_MM_sparse(std::complex<double>*  A_val_array,
					 MKL_INT* A_idx_row_array,
					 MKL_INT* A_idx_col_array,
					 MKL_INT  A_nnz,
					 std::complex<double>*  B_val_array,
					 MKL_INT  mat_dim,
					 std::complex<double>*  sol_val_array){
	
	bool print_content = false;

	_MKL_DSS_HANDLE_t dss_handle;
	_INTEGER_t        error;
	_CHARACTER_t statIn[] = "determinant";
	_DOUBLE_PRECISION_t statOut[5];

	MKL_INT  nRows     = mat_dim;
	MKL_INT  nCols     = mat_dim;
	MKL_INT  nNonZeros = A_nnz;
	MKL_INT  nRhs      = 1;
	MKL_INT* rowIndex  = A_idx_row_array;
	MKL_INT* columns   = A_idx_col_array;
	std::complex<double>*  values    = A_val_array;
	std::complex<double>*  rhs       = B_val_array;
	std::complex<double>*  solValues = sol_val_array;//[nRows];

	MKL_INT dss_opt  = MKL_DSS_MSG_LVL_WARNING + MKL_DSS_TERM_LVL_ERROR + MKL_DSS_ZERO_BASED_INDEXING;
	MKL_INT dss_sym  = MKL_DSS_NON_SYMMETRIC_COMPLEX;
	MKL_INT dss_type = MKL_DSS_POSITIVE_DEFINITE;

	/* Initialize DSS solver */
	error = dss_create(dss_handle,
					   dss_opt);
	if ( error != MKL_DSS_SUCCESS ){raise_error_DSS_MKL(error, "dss_handle");}
	if(print_content){std::cout << "Finished dss_create" << std::endl;}

	error = dss_define_structure(dss_handle,
								 dss_sym,
								 rowIndex,
								 nRows,
								 nCols,
								 columns,
								 nNonZeros);
	if ( error != MKL_DSS_SUCCESS ){raise_error_DSS_MKL(error, "dss_define_structure");}
	if(print_content){std::cout << "Finished dss_define_structure" << std::endl;}

	error = dss_reorder(dss_handle,
						dss_opt,
						0);
	if ( error != MKL_DSS_SUCCESS ){raise_error_DSS_MKL(error, "dss_reorder");}
	if(print_content){std::cout << "Finished dss_reorder" << std::endl;}

	error = dss_factor_complex(dss_handle,
							   dss_type,
							   values );
	if ( error != MKL_DSS_SUCCESS ){raise_error_DSS_MKL(error, "dss_factor_real");}
	if(print_content){std::cout << "Finished dss_factor_real" << std::endl;}

	error = dss_solve_complex(dss_handle,
						   	  dss_opt,
						   	  rhs,
						   	  nRhs,
						   	  solValues );
	if ( error != MKL_DSS_SUCCESS ){raise_error_DSS_MKL(error, "dss_solve_real");}
	if(print_content){std::cout << "Finished dss_solve_real" << std::endl;}

	if ( nRows < A_nnz ) {
		error = dss_statistics(dss_handle,
							   dss_opt,
							   statIn,
							   statOut);
		if ( error != MKL_DSS_SUCCESS ){raise_error_DSS_MKL(error, "dss_statistics");}
		
		if(print_content){
			printf(" - Determinant power is %g \n", statOut[0]);
			printf(" - Determinant base is %g \n", statOut[1]);
			printf(" - Determinant is %g \n", (pow(10.0,statOut[0]))*statOut[1]);
		}
	}
	if(print_content){std::cout << "Finished dss_statistics" << std::endl;}

	error = dss_delete(dss_handle,
					   dss_opt);
	if ( error != MKL_DSS_SUCCESS ){raise_error_DSS_MKL(error, "dss_delete");}
	if(print_content){std::cout << "Finished dss_delete" << std::endl;}

	if(print_content){
		printf(" Solution array: ");
		for(MKL_INT i = 0; i< nCols; i++)
			printf(" %g %g", solValues[i].real(), solValues[i].imag() );
		printf("\n");
	}
}

void simple_transpose_matrix_routine(double* mat_array, int mat_dim){
	double* mat_array_temp = new double [mat_dim * mat_dim];

	/* Copy to temp array */
	for (int row=0; row<mat_dim; row++){
		for (int col=0; col<mat_dim; col++){
			mat_array_temp[col*mat_dim + row] = mat_array[row*mat_dim + col];
		}
	}
	/* Copy from temp array */
	std::copy(mat_array_temp, mat_array_temp + mat_dim*mat_dim, mat_array);

	/* Delete temporary array */
	delete [] mat_array_temp;
}

void square_dense_to_sparse_COO_format_converter(int      mat_dim,
												 double*  mat_dense_array,
												 double** mat_sparse_val_array,
												 int**    mat_sparse_row_array,
												 int**    mat_sparse_col_array,
												 int&     mat_sparse_dim){

	/* Temporary arrays to hold indices of non-zero elements */
	int* temp_dense_idx_row_array = new int [mat_dim*mat_dim];
	int* temp_dense_idx_col_array = new int [mat_dim*mat_dim];
	
	/* Number of non-zero elements ("NNZ" - standard abbreviation) */
	int nnz = 0;
	
	/* Determine non-zero elements in dense-matrix and store indices to temporay array */
	for (int i=0; i<mat_dim; i++){
		for (int j=0; j<mat_dim; j++){
			if (mat_dense_array[i*mat_dim + j] != 0){
				temp_dense_idx_row_array[nnz] = i; // row-index
				temp_dense_idx_col_array[nnz] = j; // col-index

				nnz += 1;
			}
		}
	}
	
	/* Allocate minimally-sized sparse-array */
	*mat_sparse_val_array = new double [nnz];
	*mat_sparse_row_array = new int    [nnz];
	*mat_sparse_col_array = new int    [nnz];

	int row_idx = 0;
	int col_idx = 0;
	for (int i=0; i<nnz; i++){
		row_idx = temp_dense_idx_row_array[i];
		col_idx = temp_dense_idx_col_array[i];

		(*mat_sparse_val_array)[i]  = mat_dense_array[row_idx*mat_dim + col_idx];

		(*mat_sparse_row_array)[i] = row_idx;
		(*mat_sparse_col_array)[i] = col_idx;
	}
	
	mat_sparse_dim = nnz;

	delete [] temp_dense_idx_row_array;
	delete [] temp_dense_idx_col_array;
}

void square_sparse_COO_to_dense_format_converter(int      mat_dim,
												 double** mat_dense_array,
												 double*  mat_sparse_val_array,
												 int*     mat_sparse_row_array,
												 int*     mat_sparse_col_array,
												 int      mat_sparse_dim){
	
	*mat_dense_array = new double [mat_dim*mat_dim];

	/* Fill dense matrix with zeroes */
	for (int i=0; i<mat_dim; i++){
		mat_dense_array[i] = 0;
	}

	int row_idx = 0;
	int col_idx = 0;
	/* Write sparse non-zero elements into dense matrix */
	for (int i=0; i<mat_sparse_dim; i++){
		row_idx = mat_sparse_row_array[i];
		col_idx = mat_sparse_col_array[i];

		(*mat_dense_array)[row_idx*mat_dim + col_idx] = mat_sparse_val_array[i];
	}
}

void coo_to_csr_format_converter(int* idx_row_array_coo,
								 int* idx_row_array_csr,
								 int  mat_sparse_dim,
								 int  mat_dense_dim){
	for (int i=0; i<mat_dense_dim+1; i++){
		idx_row_array_csr[i] = 0;
	}

	for (int i=1; i<mat_dense_dim+1; i++){

		int num_nnz_cols_in_row = 0;
		for (int j=idx_row_array_csr[i-1]; j<mat_sparse_dim; j++){
			if (idx_row_array_coo[j]>i-1){
				break;
			}
			else{
				num_nnz_cols_in_row += 1;
			}
		}

		idx_row_array_csr[i] = num_nnz_cols_in_row + idx_row_array_csr[i-1];
	}
}

void coo_col_major_to_coo_row_major_converter(double** mat_val_array_coo,
											  int**    mat_row_array_coo,
											  int**    mat_col_array_coo,
											  int      mat_sparse_dim,
											  int      mat_dense_dim){
	
	/* Allocate row-major (RM) arrays (these will replace input arrays) */
	double* mat_val_array_coo_RM_temp = new double [mat_sparse_dim];
	int* 	mat_row_array_coo_RM_temp = new int    [mat_sparse_dim];
	int* 	mat_col_array_coo_RM_temp = new int    [mat_sparse_dim];

	/* Loop through RM row-array */
	int nnz_counter = 0;
	for (int row_RM=0; row_RM<mat_dense_dim; row_RM++){
		for (int col_RM=0; col_RM<mat_dense_dim; col_RM++){
			
			/* Random order loop */
			for (int nnz_idx=0; nnz_idx<mat_sparse_dim; nnz_idx++){
				int row_CM = (*mat_row_array_coo)[nnz_idx];
				int col_CM = (*mat_col_array_coo)[nnz_idx];

				/* If nnz element corresponds to order of two outer loops, append element */
				if (row_RM==row_CM and col_RM==col_CM){
					/* Append element */
					mat_row_array_coo_RM_temp[nnz_counter] = row_RM;
					mat_col_array_coo_RM_temp[nnz_counter] = col_RM;
					mat_val_array_coo_RM_temp[nnz_counter] = (*mat_val_array_coo)[nnz_idx];

					/* Increment non-zero counter */
					nnz_counter += 1;

					/* Break innermost loop */
					break;
				}
			}
		}
	}

	/* Delete old column-major format arrays */
	delete [] *mat_val_array_coo;
	delete [] *mat_row_array_coo;
	delete [] *mat_col_array_coo;

	/* Let input array pointers point to temp-arrays */
	*mat_val_array_coo = mat_val_array_coo_RM_temp;
	*mat_row_array_coo = mat_row_array_coo_RM_temp;
	*mat_col_array_coo = mat_col_array_coo_RM_temp;
}

void increase_sparse_array_size(double** sparse_array, int array_length, int sparse_step_length){
	/*Allocate an array to hold sparse_array values with increased size */
	double* sparse_array_temp = new double [array_length + sparse_step_length];

	/* Copy elements from sparse_array to sparse_array_temp */
	std::copy(*sparse_array, *sparse_array + array_length, sparse_array_temp);

	/* Free memory pointed to by sparse_array */
	delete [] *sparse_array;

	/* Let sparse_array point to newly created, expanded array sparse_array_temp */
	*sparse_array = sparse_array_temp;
}
void increase_sparse_array_size(int** sparse_array, int array_length, int sparse_step_length){
	/*Allocate an array to hold sparse_array values with increased size */
	int* sparse_array_temp = new int [array_length + sparse_step_length];

	/* Copy elements from sparse_array to sparse_array_temp */
	std::copy(*sparse_array, *sparse_array + array_length, sparse_array_temp);

	/* Free memory pointed to by sparse_array */
	delete [] *sparse_array;

	/* Let sparse_array point to newly created, expanded array sparse_array_temp */
	*sparse_array = sparse_array_temp;
}

void reduce_sparse_array_size(double** sparse_array, int array_length, int sparse_dim){
	/*Allocate an array to hold sparse_array values with minimal size */
	double* sparse_array_temp = new double [sparse_dim];

	/* Copy elements from sparse_array to sparse_array_temp */
	std::copy(*sparse_array, *sparse_array + sparse_dim, sparse_array_temp);

	/* Free memory pointed to by sparse_array */
	delete [] *sparse_array;
	
	/* Let sparse_array point to newly created, minimal-size array sparse_array_temp */
	*sparse_array = sparse_array_temp;
}
void reduce_sparse_array_size(int** sparse_array, int array_length, int sparse_dim){
	/*Allocate an array to hold sparse_array values with minimal size */
	int* sparse_array_temp = new int [sparse_dim];

	/* Copy elements from sparse_array to sparse_array_temp */
	std::copy(*sparse_array, *sparse_array + sparse_dim, sparse_array_temp);

	/* Free memory pointed to by sparse_array */
	delete [] *sparse_array;
	
	/* Let sparse_array point to newly created, minimal-size array sparse_array_temp */
	*sparse_array = sparse_array_temp;
}
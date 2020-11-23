
#include "lippmann_schwinger_solver.h"

/* A -> A*B, where B is a scalar and A a vector */
void cdot_VS_colMaj(std::complex<float>* A, std::complex<float>* B, int N){
	cblas_cscal(N, B, A, N);
}
void cdot_VS_colMaj(std::complex<double>* A, std::complex<double>* B, int N){
	cblas_zscal(N, B, A, N);
}
void solve_MM_lin_eq(std::complex<float>* A, std::complex<float>* B, int N){
	char trans = 'N';
	long long int ipiv [N];
	
	LAPACKE_cgetrf(LAPACK_ROW_MAJOR, N, N, A, N, ipiv);
    LAPACKE_cgetrs(LAPACK_ROW_MAJOR, trans, N, N, A, N, ipiv, B, N);
}
void solve_MM_lin_eq(std::complex<double>* A, std::complex<double>* B, int N){
	char trans = 'N';
	long long int ipiv [N];
	
	LAPACKE_zgetrf(LAPACK_ROW_MAJOR, N, N, A, N, ipiv);
    LAPACKE_zgetrs(LAPACK_ROW_MAJOR, trans, N, N, A, N, ipiv, B, N);
}
void M_invert(std::complex<float>* A, int N){
	long long int ipiv [N];

	//int info = 0;
	//int lwork = N*N;
	//std::complex<float> work [lwork];
	
	LAPACKE_cgetrf(LAPACK_ROW_MAJOR, N, N, A, N, ipiv);
    LAPACKE_cgetri(LAPACK_ROW_MAJOR, N, A, N, ipiv);
}
void M_invert(std::complex<double>* A, int N){
	long long int ipiv [N];

	//int info = 0;
	//int lwork = N*N;
	//std::complex<double> work [lwork];
	
	LAPACKE_zgetrf(LAPACK_ROW_MAJOR, N, N, A, N, ipiv);
    LAPACKE_zgetri(LAPACK_ROW_MAJOR, N, A, N, ipiv);
}
/* Dot product for C := alpha*A*B + beta*C for dimensions C(N,M), A(N,K), and B(K,M) (row-major) */
void cdot_MM(std::complex<float> *A, std::complex<float> *B, std::complex<float> *C, int N, int K, int M){
	float beta = 0;
	float alpha = 1;
	cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, M, K, &alpha, A, K, B, M, &beta, C, M);
}
void cdot_MM(std::complex<double> *A, std::complex<double> *B, std::complex<double> *C, int N, int K, int M){
	double beta = 0;
	double alpha = 1;
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, M, K, &alpha, A, K, B, M, &beta, C, M);
}


void make_denominator_array(cfloatType* D_array, int Nk, double* k_array, double* wk_array, double E, double M){

	const cfloatType J{ 0.0, 1.0 };

	bool E_positive = (E>=0);

	/* Reset D */
	for (int i=0; i<Nk+1; i++){
		D_array[i] = 0;
	}
	
	floatType norm_fac = 1;
	//floatType norm_fac = 2/M_PI;

	double temp = 0;
	for (int i=0; i<Nk; i++){
		//temp = norm_fac*wk_array[i] / (E - k_array[i]*k_array[i]);
		temp = norm_fac*wk_array[i] / (k_array[i]*k_array[i]/MN - E);

		/* LS integral summation */
		D_array[i]  = temp * k_array[i] * k_array[i];

		if (E_positive){
			/* Principle value prescription term */
			D_array[Nk] -= temp;
		}
	}

	if (E_positive){
		floatType q = sqrt(E*M);
		D_array[Nk] *= q*q;

		/* Complex integration term. Only required in on-shell calculations */
		floatType mu = 0.5*M;
		D_array[Nk] += J*q*M_PI*mu;
	}
}

void make_wave_matrix(cfloatType* F_array, cfloatType* D_array, int Nk1, bool coupled){

	for (int i=0; i<Nk1; i++){
		cfloatType factor = D_array[i];
		
		if (coupled == false){
			/* Calculate V*D */
			cdot_VS_colMaj(F_array + i, &factor, Nk1);

			/* Add identity matrix */
			F_array[i*Nk1 + i] += 1;
		}
		else{
			/* Calculate V*D */
			cdot_VS_colMaj(F_array + i	    , &factor, 2*Nk1);
			cdot_VS_colMaj(F_array + i + Nk1, &factor, 2*Nk1);

			/* Add identity matrix */
			F_array[i*2*Nk1 + i] 		     += 1;
			F_array[(i+Nk1)*2*Nk1 + i + Nk1] += 1;
		}
	}
}

void calculate_t_element(double* V_prestored_array,
						   cfloatType* T_array,
						   bool coupled,
						   double E, double M,
						   int Nk, double* k_array, double* wk_array,
						   int idx_row, int idx_col){

	/* Number of quadrature points plus on-shell momentum (x2 for coupled matrices) */
	int Nk1 = Nk+1;
	int Nk2 = 2*Nk1;

	/* Matrix dimension */
	int mat_dim = Nk1;

	/* Detemine if this is a coupled channel */
	if (coupled){
		mat_dim *= 2;
	}
	int mat_dim_sq = mat_dim * mat_dim;

	floatType* 	V_array = V_prestored_array;			// Potential matrix
	cfloatType* D_array = new cfloatType [Nk1];			// Denominator vector
	cfloatType* F_array = new cfloatType [mat_dim_sq];	// Wave matrix
	//cfloatType* T_array = new cfloatType [mat_dim_sq];	// T-matrix

	/* Set resolvent array pointer to pre-calculated arrays */
	make_denominator_array(D_array, Nk, k_array, wk_array, E, M);
	
	/* Copy potentials into wave- and T-matrices
	 * (required since to BLAS-routines rewrites matrices). */
	for (int i=0; i<mat_dim_sq; i++){
		F_array[i] = V_array[i];
		T_array[i] = V_array[i];
	}
	
	/* Calculate wave-matrix */
	make_wave_matrix(F_array, D_array, Nk1, coupled);
	
	/* Solve LS */
	solve_MM_lin_eq(F_array, T_array, mat_dim);

	//M_invert(F_array, mat_dim);
	//cfloatType* V_array_complex_cast = new cfloatType [mat_dim_sq];
	//for (int i=0; i<mat_dim_sq; i++){
	//	V_array_complex_cast[i] = V_array[i];
	//}
	//cdot_MM(F_array, V_array_complex_cast, T_array, Nk1, Nk1, Nk1);
	//delete [] V_array_complex_cast;


/*
	if (coupled){
		for (int i=0; i<Nk; i++){
			for (int j=0; j<Nk; j++){
				std::cout << T_array[i*Nk2 + j] << std::endl; 
				std::cout << T_array[(i+Nk1)*Nk2 + j] << std::endl; 
				std::cout << T_array[i*Nk2 + j+Nk1] << std::endl; 
				std::cout << T_array[(i+Nk1)*Nk2 + j+Nk1] << std::endl; 
			}
		}
	}
*/

	//for (int i=0; i<mat_dim_sq; i++){
	//	t_array[i] = T_array[i].real();
	//}
	
	/* Delete temporary wave matrix */
	//delete [] V_array;
	delete [] D_array;
	delete [] F_array;
	//delete [] T_array;
}
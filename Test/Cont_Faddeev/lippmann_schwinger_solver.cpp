
#include "lippmann_schwinger_solver.h"

/* A -> A*B, where B is a scalar and A a vector */
void cdot_VS_colMaj(std::complex<float> *A, std::complex<float> *B, int N){
	cblas_cscal(N, B, A, N);
}
void cdot_VS_colMaj(std::complex<double> *A, std::complex<double> *B, int N){
	cblas_zscal(N, B, A, N);
}

void make_denominator_array(cfloatType *D){
	
	/* Reset D */
	for (int i=0; i<N1; i++){
		D[i] = 0;
	}
	
	floatType temp;
	for (int i=0; i<N; i++){
		temp = (*wList)[i] * M / ((*qList)[i]*(*qList)[i] - (*qList)[N]*(*qList)[N]);
		D[i]  = temp * (*qList)[i] * (*qList)[i];
		D[N] -= temp;
	}
	
	D[N] *= (*qList)[N]*(*qList)[N];
	D[N] += divTwo*J*(*qList)[N]*pi*M;
}

void make_wave_matrix(cfloatType *F, cfloatType *D, bool isCoupled){
	
	for (int i=0; i<N1; i++){
		cfloatType factor = D[i];
		
		if (isCoupled == false){
			cdot_VS_colMaj(F+i, &factor, N1);
			F[N1*i+i] += 1;
		}
		else{
			cdot_VS_colMaj(F+i, &factor, 2*N1);
			cdot_VS_colMaj(F+i+N1, &factor, 2*N1);
			F[2*N1*i+i] 		+= 1;
			F[2*N1*(i+N1)+i+N1] += 1;
		}
	}
}

void solve_for_T_element(){
    /* On-shell momentum */
	floatType q = (*expQcom)[I];
	(*qList)[N] = q;
	
	/* Update potential matrix */
	potentialMatrix->updateMatrix(V, L, S, J, q);
	
	/* Set resolvent array pointer to pre-calculated arrays */
	makeDenomVec(D);
	
	cfloatType *F = new cfloatType [N2];
	
	/* Copy potentials into wave- and T-matrices
	 * (needed due to BLAS-routines rewriting
	 * matrices). */
	for (int i=0; i<N2; i++){
		F[i] = V[i];
		T[i] = V[i];
	}
	
	/* Calculate wave-matrices */
	waveMatrixCBLAS(F, D, false);
	
	/* Solve LS */
	solveMM(F, T, N1);
	
	/* Append on-shell elements */
	Telements[I] = T[N2 - 1];
	
	/* Delete temporary wave matrix */
	delete [] F;
}

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

/*double extract_potential_element_from_array(int L, int Lp, int J, int S, bool coupled, double* V_array){
    
	if (J<abs(S-L) or J>S+L or J<abs(S-Lp) or J>S+Lp){
		std::cout << L <<" "<< Lp <<" "<< J <<" "<< S << std::endl;
		raise_error("Encountered unphysical state in LS-solver");
	}

    double potential_element = NAN;

    if (coupled){
        if (L==Lp and L<J){         // --
            potential_element = V_array[2];
        }
        else if (L==Lp and L>J){    // ++
            potential_element = V_array[5];
        }
        else if (L<Lp){             // +-
            potential_element = V_array[3];
        }
        else{                       // -+
            potential_element = V_array[4];
        }
    }
    else{
        if (J==0 and L!=J){         // 3P0 (++)
            potential_element = V_array[5];
        }
        else if (S==0){             // S=0
            potential_element = V_array[0];
        }
        else{                       // S=1
            potential_element = V_array[1];
        }
    }

    return potential_element;
}*/

void make_denominator_array(cfloatType* D_array, int Nk, double* k_array, double* wk_array, double E, double M){

	const cfloatType J{ 0.0, 1.0 };

	bool E_positive = (E>=0);

	/* Reset D */
	for (int i=0; i<Nk+1; i++){
		D_array[i] = 0;
	}
	
	double norm_fac = 1;
	//double norm_fac = 2/M_PI;

	double temp = 0;
	for (int i=0; i<Nk; i++){
		temp = norm_fac*wk_array[i] / (k_array[i]*k_array[i]/M - E);
		/* LS integral summation */
		D_array[i]  = temp * k_array[i] * k_array[i];

		if (E_positive){
			/* Principle value prescription term */
			D_array[Nk] -= temp;
		}
	}

	if (E_positive){
		double q = sqrt(E*M);
		D_array[Nk] *= q*q;

		/* Complex integration term. Only required in on-shell calculations */
		D_array[Nk] += 0.5*J*q*pi*M;
	}
	else{
		D_array[Nk] = 0;
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

/* Calculates half off-shell matrix t(q,p;E) where E=q^2/M (i.e. q is on-shell).
 * Returns matrix element given by idx_row and idx_col */
double calculate_t_element(double* V_prestored_array,
						   double* t_array,
						   int L, int L_p, int S, int J, int T,
						   double E, double M,
						   int Nk, double* k_array, double* wk_array,
						   int idx_row, int idx_col,
						   potential_model* pot_ptr_nn,
						   potential_model* pot_ptr_np){

	/* Number of quadrature points plus on-shell momentum (x2 for coupled matrices) */
	int Nk1 = Nk+1;
	int Nk2 = 2*Nk1;

	/* Matrix dimension */
	int mat_dim = Nk1;

	/* Detemine if this is a coupled channel */
	bool coupled = false;
	if (L!=L_p or (L==L_p and L!=J and J!=0)){
		coupled  = true;
		mat_dim *= 2;
	}
	int mat_dim_sq = mat_dim * mat_dim;

	double* 	V_array = new double 	 [mat_dim_sq];								// Potential matrix
	cfloatType* D_array = new cfloatType [Nk1];		// Denominator vector
	cfloatType* F_array = new cfloatType [mat_dim_sq];	// Wave matrix
	cfloatType* T_array = new cfloatType [mat_dim_sq];	// T-matrix

	
	for (int i=0; i<Nk; i++){
		/* Copy matrix quadrature gridpoints */
		for (int j=0; j<Nk; j++){
			if (coupled){
				V_array[ i	   *Nk2 +  j] 	   = V_prestored_array[ i	 *2*Nk +  j];
				V_array[(i+Nk1)*Nk2 +  j] 	   = V_prestored_array[(i+Nk)*2*Nk +  j];
				V_array[ i	   *Nk2 + (j+Nk1)] = V_prestored_array[ i	 *2*Nk + (j+Nk)];
				V_array[(i+Nk1)*Nk2 + (j+Nk1)] = V_prestored_array[(i+Nk)*2*Nk + (j+Nk)];
			}
			else{
				V_array[i*Nk1 + j] = V_prestored_array[i*Nk + j];
			}
		}

		/* Copy elements for end row and column */
		if (coupled){
			/* Set end-column <ki|v|p> */
			V_array[ i	   *Nk2 +  Nk] 		= V_prestored_array[ i	  *2*Nk +  idx_col];
			V_array[(i+Nk1)*Nk2 +  Nk] 		= V_prestored_array[(i+Nk)*2*Nk +  idx_col];
			V_array[ i	   *Nk2 + (Nk+Nk1)] = V_prestored_array[ i	  *2*Nk + (idx_col+Nk)];
			V_array[(i+Nk1)*Nk2 + (Nk+Nk1)] = V_prestored_array[(i+Nk)*2*Nk + (idx_col+Nk)];
			/* Set end-row <p'|v|ki> */
			V_array[ Nk		*Nk2 +  i] 		= V_prestored_array[ idx_row	*2*Nk +  i];
			V_array[(Nk+Nk1)*Nk2 +  i] 		= V_prestored_array[(idx_row+Nk)*2*Nk +  i];
			V_array[ Nk		*Nk2 + (i+Nk1)] = V_prestored_array[ idx_row	*2*Nk + (i+Nk)];
			V_array[(Nk+Nk1)*Nk2 + (i+Nk1)] = V_prestored_array[(idx_row+Nk)*2*Nk + (i+Nk)];
		}
		else{
			/* Set end-column <ki|v|p> */
			V_array[i*mat_dim + Nk] = V_prestored_array[i*Nk + idx_col];
			/* Set end-row <p'|v|ki> */
			V_array[Nk*mat_dim + i] = V_prestored_array[idx_row*Nk + i];
		}
	}
	if (coupled){
		V_array[ Nk		*Nk2 +  Nk] 	 = V_prestored_array[ idx_row	 *2*Nk +  idx_col];
		V_array[(Nk+Nk1)*Nk2 +  Nk] 	 = V_prestored_array[(idx_row+Nk)*2*Nk +  idx_col];
		V_array[ Nk		*Nk2 + (Nk+Nk1)] = V_prestored_array[ idx_row	 *2*Nk + (idx_col+Nk)];
		V_array[(Nk+Nk1)*Nk2 + (Nk+Nk1)] = V_prestored_array[(idx_row+Nk)*2*Nk + (idx_col+Nk)];
	}
	else{
		/* Set last element <p'|v|p> */
		V_array[mat_dim_sq - 1] = V_prestored_array[idx_row*Nk + idx_col];
	}


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
	
	cfloatType T_element = NAN;
	if (coupled){
		if (L>L_p){		 // upper-right block of coupled T-matrix
			T_element = T_array[idx_row*Nk2 + idx_col + Nk1]; 
		}
		else if (L<L_p){ // lower-left block of coupled T-matrix
			T_element = T_array[(idx_row + Nk1)*Nk2 + idx_col];
		}
		else if (L<J){	 // upper-left block of coupled T-matrix
			T_element = T_array[idx_row*Nk1 + idx_col];
		}
		else{			 // lower-right block of coupled T-matrix
			T_element = T_array[(idx_row + Nk1)*Nk2 + idx_col + Nk1];
		}
	}
	else{
		T_element = T_array[idx_row*Nk1 + idx_col];
	}

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

	for (int i=0; i<mat_dim_sq; i++){
		t_array[i] = T_array[i].real();
	}

	double T_element_real = T_element.real();
	
	/* Delete temporary wave matrix */
	delete [] V_array;
	delete [] D_array;
	delete [] F_array;
	delete [] T_array;

	return T_element_real;
}
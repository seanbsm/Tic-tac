
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

double extract_potential_element_from_array(int L, int Lp, int J, int S, bool coupled, double* V_array){
    if (J<abs(S-L) or J>S+L){exit(1);}
    if (J<abs(S-Lp) or J>S+Lp){exit(1);}

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
}

void make_denominator_array(cfloatType* D_array, int Nk, double* k_array, double* wk_array, double q, double M){

	/* Reset D */
	for (int i=0; i<Nk+1; i++){
		D_array[i] = 0;
	}
	
	double temp = 0;
	for (int i=0; i<Nk; i++){
		temp = M*wk_array[i] / (k_array[i]*k_array[i] - q*q);
		/* LS integral summation */
		D_array[i]  = temp * k_array[i] * k_array[i];
		/* Principle value prescription term */
		D_array[Nk] -= temp;
	}
	D_array[Nk] *= q*q;

	/* Complex integration term. Only required in on-shell calculations */
	D_array[Nk] += 0.5*J*q*pi*M;
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
			cdot_VS_colMaj(F_array + i		 , &factor, 2*Nk1);
			cdot_VS_colMaj(F_array + i + Nk1,  &factor, 2*Nk1);

			/* Add identity matrix */
			F_array[i*2*Nk1 + i] 			 += 1;
			F_array[(i+Nk1)*2*Nk1 + i + Nk1] += 1;
		}
	}
}

/* Calculates half off-shell matrix t(q,p;E) where E=q^2/M (i.e. q is on-shell).
 * Returns matrix element given by idx_row and idx_col */
double calculate_t_element(int L, int L_p, int S, int J, int T,
						   double E, double p, double M,
						   int Nk, double* k_array, double* wk_array,
						   int idx_row, int idx_col,
						   potential_model* pot_ptr_nn,
						   potential_model* pot_ptr_np){

	double q = sqrt(E*M);

	/* Number of quadrature points plus on-shell momentum (x2 for coupled matrices) */
	int Nk1 = Nk + 1;
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
	
	/* Set resolvent array pointer to pre-calculated arrays */
	make_denominator_array(D_array, Nk, k_array, wk_array, q, M);
	
	/* Construct 2N potential matrix <k|v|k_p> */
	double V_IS_elements [6];	// Isoscalar (IS)
	double V_nn_elements [6];	// neutron-neutron (nn)
	double V_np_elements [6];	// neutron-proton (np)
	double k=0, k_in=0, k_p=0, k_out=0;
	for (int idx_k=0; idx_k<Nk1; idx_k++){
        if (idx_k<Nk){k = k_array[idx_k];}	// We set to quadrature momentum
		else		 {k = p;}				// We set to off-shell momentum p
        k_in = k;//*hbarc;
                    
        for (int idx_k_p=0; idx_k_p<Nk1; idx_k_p++){
			if (idx_k_p<Nk){k_p = k_array[idx_k_p];}	// We set to quadrature momentum
			else		   {k_p = q;} 					// We set to on-shell momentum q
            k_out = k_p;//*hbarc;

			/* We create an isoscalar potential */
			if (T==1){ // Interaction can be either nn or np
                pot_ptr_nn->V(k_in, k_out, coupled, S, J, T, V_nn_elements);
                pot_ptr_np->V(k_in, k_out, coupled, S, J, T, V_np_elements);

				for (int idx_element=0; idx_element<6; idx_element++){
					V_IS_elements[idx_element] = (1./3)*V_np_elements[idx_element] + (2./3)*V_nn_elements[idx_element];
				}
            }
            else{ 	   // Interaction must be np
                pot_ptr_np->V(k_in, k_out, coupled, S, J, T, V_IS_elements);
            }

			/* Write element to potential matrix V_array */
			if (coupled){
				V_array[idx_k*Nk2 + idx_k_p]			   = extract_potential_element_from_array(J-1, J-1, J, S, coupled, V_IS_elements);
				V_array[idx_k*Nk2 + idx_k_p + Nk1]		   = extract_potential_element_from_array(J-1, J+1, J, S, coupled, V_IS_elements);
				V_array[(idx_k + Nk1)*Nk2 + idx_k_p]	   = extract_potential_element_from_array(J+1, J-1, J, S, coupled, V_IS_elements);
				V_array[(idx_k + Nk1)*Nk2 + idx_k_p + Nk1] = extract_potential_element_from_array(J+1, J+1, J, S, coupled, V_IS_elements);
			}
			else{
				V_array[idx_k*Nk1 + idx_k_p] = extract_potential_element_from_array(L, L_p, J, S, coupled, V_IS_elements);
			}
		}
	}
	
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
			T_element = T_array[idx_row*Nk2 + idx_col];
		}
		else{			 // lower-right block of coupled T-matrix
			T_element = T_array[(idx_row + Nk1)*Nk2 + idx_col + Nk1];
		}
	}
	else{
		T_element = T_array[idx_row*Nk1 + idx_col];
	}

	double T_element_real = T_element.real();

	std::cout << T_element << std::endl;

	/* Delete temporary wave matrix */
	delete [] V_array;
	delete [] D_array;
	delete [] F_array;
	delete [] T_array;

	return T_element_real;
}
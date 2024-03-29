
#include "faddeev_iterator.h"


/* Finds the eigenvalues of a real,
 * general N-by-N matrix A using LAPACK: geev.
 * "wr" will be filled with the real components of
 * the eigenvalues, in ascending order. We work
 * with row-major matrices, as is usual
 * in C & C++ */
void find_eigenvalues(double* A, double* wr, double* vr, int N){
    
    char jobvl = 'N'; // left eigenvectors (u in u*A=u*w) of A are not computed.
    char jobvr = 'V'; // right eigenvectors (v in A*v=w*v) of A are computed. 

    /* Imaginary part of eigenvalues.
     * This variable is not used anywhere
     * and is merely required by the routine */
    double  wi [N];

    double* vl = NULL; // This will not be referenced since jobvl='N'

    int lda  = N;   // leading dimension of A
    int ldvr = N;   // leading dimension of vr
    int ldvl = N;   // leading dimension of vl

    LAPACKE_dgeev(LAPACK_ROW_MAJOR, jobvl, jobvr, N, A, lda, wr, wi, vl, ldvl, vr, ldvr);
}


/* Copied from https://stackoverflow.com/questions/25313816/gram-schmidt-function-not-working-c/25320243
 * The function uses the modified Gram-Scmidt orgonalization routine to read a matrix of vectors (state_matrix)
 * to create an orgothogonal basis from them (state_basis).
 * Each row corresponds to a state.  */
void modified_gram_schmidt(double* state_matrix,
                           double* state_basis,
                           int num_states,
                           int num_state_elements){
    
    /* Define easier-to-read internal notation */
    int rows = num_states;
    int cols = num_state_elements;

    double v [rows*cols];

    /* Copy input state_matrix into v */
    for (int i=0; i<rows; i++){   // state loop
        for (int j=0; j<cols; j++){ // element loop
            v[i*cols + j] = state_matrix[i*cols + j];
        }
    }

    /* Input state normalisation */
    double norm;
    double dot_product;

    /* Loop over states */
    for (int i=0; i<rows; i++){        // state loop

        /* Calculate normalisation of vi */
        norm = 0;
        for (int n=0; n<cols; n++){    // element loop
            norm += v[i*cols + n] * v[i*cols + n];
        }
        norm = sqrt(norm);

        /* Set orthogonal basis vector i to vi/|vi| */
        for (int j=0; j<cols; j++){       // element loop
            state_basis[i*cols + j] = v[i*cols + j] / norm;
        }

        /* Set basis states k orthogonal to input state i */
        for (int k=i+1;  k<rows; k++){ // state loop

            /* Calculate dot-product of state i in state_basis with vk */
            dot_product = 0;
            for (int n=0; n<cols; n++){
                dot_product += state_basis[i*cols + n] * v[k*cols + n];
            }
            
            /* Modify state vk */
            for (int j=0; j<cols; j++){ // element loop
                v[k*cols + j] = v[k*cols + j] - dot_product * state_basis[i*cols + j];
            }
        }
    }
}

void brute_force_lanczos_for_faddeev(double* states_array,
                                     double* physical_state_array,
                                     int& physical_state_idx,
                                     int num_states,
                                     int num_state_elements){
    
    /* Orthogonal state basis */
    double ort_states_array [num_states*num_state_elements];

    /* Call modified Gram-Schmidt routine to fill ort_states_array */
    modified_gram_schmidt(states_array,
                          ort_states_array,
                          num_states,
                          num_state_elements);

    /* Calculate a and b state-bases transformation matrices */
    double a_array [num_states*num_states];   // Transforms ort_state to state
    double b_array [num_states*num_states];   // Transforms state to ort_state
    double inner_product_for_a = NAN;
    double inner_product_for_b = NAN;
    for (int i=0; i<num_states; i++){
        /* a and b are lower-triangular due to Gram-Schmidt orthogonalization */
        for (int j=0; j<i+1; j++){
            inner_product_for_a = 0;
            inner_product_for_b = 0;

            /* Inner-products of states */
            for (int k=0; k<num_state_elements; k++){
                inner_product_for_a += ort_states_array[j*num_state_elements + k] * states_array[i*num_state_elements + k];
                inner_product_for_b += states_array[j*num_state_elements + k] * ort_states_array[i*num_state_elements + k];
            }
            
            a_array[i*num_states + j] = inner_product_for_a;
            b_array[i*num_states + j] = inner_product_for_b;
        }
    }

    /* Calculate M-matrix from a and b */
    double M_array [num_states*num_states];
    for (int i=0; i<num_states; i++){
        for (int n=0; n<num_states; n++){
            
            M_array[i*num_states + n] = 0;
            for (int j=0; j<i+1; j++){
                M_array[i*num_states + n] += b_array[i*num_states + j] * a_array[(j+1)*num_states + n];
            }
        }
    }

    /* Diagonalize M */
    double M_eigenvalues  [num_states];
    double M_eigenvectors [num_states*num_states];
    find_eigenvalues(M_array, M_eigenvalues, M_eigenvectors, num_states);
    
    double closest_eigenvalue_to_one = 0;
    /* See if any an eigenvalue is close to 1,
     * and set state_value to this eigenvector */
    for (int i=0; i<num_states; i++){
        /* See if we have a physical state (must select eigenvalue closest to 1) */
        if (abs(M_eigenvalues[i]-1)<1e-14 and closest_eigenvalue_to_one<M_eigenvalues[i] and M_eigenvalues[i]<=1){
            physical_state_idx = i;
        }
    }

    /* Express physical state in orthogonal basis */
    /* Loop over linear expansion in ort_states */
    for (int j=0; j<num_states; j++){
        /* Loop over elements in given ort_state */
        for (int k=0; k<num_state_elements; k++){
            physical_state_array[k] += M_eigenvectors[physical_state_idx*num_states+j] * ort_states_array[j*num_state_elements + k];
        }
    }
}

void make_G_array(double* G_array,
                  int Np, double* p_array,
                  int Nq, double* q_array,
                  int Nx, double* x_array,
                  int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array, int* l_3N_array, int* two_j_3N_array,
                  int two_T, int two_J){
    
    /* Change in notation from my program to script in this function */
    int     Jj_dim       = Nalpha;
    int     Nx_Gtilde    = Nx;
    int     Np_3N        = Np;
    int     Nq_3N        = Nq;
    int*    L12_Jj       = l_3N_array;
    int*    l3_Jj        = L_2N_array;
    int*    J12_Jj       = J_2N_array;
    int*    two_j3_Jj    = two_j_3N_array;
    int*    S12_Jj       = S_2N_array;
    int*    T12_Jj       = T_2N_array;
    double* p_3N         = p_array;
    double* q_3N         = q_array;
    double* x_Gtilde     = x_array;
    double* Gtilde_store = G_array;
    /* The remaining function is copied from an external script and utilises auxiliary.cpp functionality
     * (I have commented out print-commands) */
    
    // determine optimized Lmax: Lmax = max(get_L)+max(get_l)
    int max_L12 = 0;
    int max_l3 = 0;
    int max_J12 = 0;

    for (int alpha = 0; alpha <= Jj_dim - 1; alpha++)
    {
        if (J12_Jj[alpha] > max_J12) max_J12 = J12_Jj[alpha];
        if (L12_Jj[alpha] > max_L12) max_L12 = L12_Jj[alpha];
        if (l3_Jj[alpha] > max_l3) max_l3 = l3_Jj[alpha];
    }

    // for F_local_matrix prestorage and F_interpolate
    int lmax = GSL_MAX_INT(max_l3, max_L12) + 3; // for C4 it is possible to couple l=lmax with THREE Y_{1}^{mu}
    //int l_interpolate_max = 2 * (lmax - 3) + 3;

    //cout << "lmax = " << lmax << ", l_interpolate_max = " << l_interpolate_max << "\n";

    int Lmax = max_L12 + max_l3;
    //int kLegendremax = 2 * max_L12;

    //int two_jmax_Clebsch = 2 * lmax;
    //int jmax_Clebsch = lmax;
    int two_jmax_SixJ = 2 * lmax; // do we need to prestore 6j??


    // for angular integration in Gtilde
    //double x_Gtilde[Nx_Gtilde];
    //double wx_Gtilde[Nx_Gtilde];

    //calc_gauss_points (x_Gtilde, wx_Gtilde, -1.0, 1.0, Nx_Gtilde);

    //MKL_INT64 g_N = (kLegendremax + 1) * (Lmax + 1) * (Lmax + 1) * Jj_dim * Jj_dim;
    //double *gtilde_array = new double[g_N];

    //MKL_INT64 Gtilde_N = Jj_dim * Jj_dim * Np_3N * Nq_3N * Nx_Gtilde;
    MKL_INT64 Atilde_N = Jj_dim * Jj_dim * (Lmax + 1);

    //double *Gtilde_store = new double[Gtilde_N];
    double *Atilde_store = new double[Atilde_N];

    for (MKL_INT64 i = 0; i <= Atilde_N - 1; i++){
        Atilde_store[i] = 0.0;
    }

    //cout << "prestore SixJ...\n";
    double *SixJ_array = new double[(two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1)];
    //cout << "SixJ_test (>0?):" << (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) << "\n";

    std::cout << "   - Pre-storing Wigner-6j coupling symbols ..." << std::endl;
    #pragma omp parallel for collapse(3)
    
    for (int two_l1 = 0; two_l1 <= two_jmax_SixJ; two_l1++){
        for (int two_l2 = 0; two_l2 <= two_jmax_SixJ; two_l2++){
            for (int two_l3 = 0; two_l3 <= two_jmax_SixJ; two_l3++){
                for (int two_l4 = 0; two_l4 <= two_jmax_SixJ; two_l4++){
                    for (int two_l5 = 0; two_l5 <= two_jmax_SixJ; two_l5++){
                        for (int two_l6 = 0; two_l6 <= two_jmax_SixJ; two_l6++){
                            SixJ_array[
                                two_l1 * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1)
                                + two_l2 * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1)
                                + two_l3 * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1)
                                + two_l4 * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1)
                                + two_l5 * (two_jmax_SixJ + 1)
                                + two_l6
                            ]
                            // implement checks for quantum numbers because of bug in gsl library
                                = gsl_sf_coupling_6j(two_l1, two_l2, two_l3, two_l4, two_l5, two_l6);
                        }
                    }
                }
            }
        }
    }

    std::cout << "   - Calculating Atilde ..." << std::endl;
    for (int alpha = 0; alpha <= Jj_dim - 1; alpha++){
        for (int alphaprime = 0; alphaprime <= Jj_dim - 1; alphaprime++){
            for (int Ltotal = 0; Ltotal <= Lmax; Ltotal++){
                Atilde_store[alpha * Jj_dim * (Lmax + 1) + alphaprime * (Lmax + 1) + Ltotal] = Atilde (alpha, alphaprime, Ltotal, Jj_dim, L12_Jj, l3_Jj, J12_Jj, two_j3_Jj, S12_Jj, T12_Jj, two_J, two_T, SixJ_array, two_jmax_SixJ);
            }
        }
    }

    /* Calculate G_array */
    //long int fullsize = Np_3N*Nq_3N*Nx_Gtilde*Jj_dim*(Jj_dim - 1);
    //long int counter = 0;
    //int frac_n, frac_o=0;
    std::cout << "   - Calculating Gtilde ..." << std::endl;
    #pragma omp parallel
    {
        #pragma omp for

        for (MKL_INT64 p_index = 0; p_index <= Np_3N - 1; p_index++){
            for (MKL_INT64 q_index = 0; q_index <= Nq_3N - 1; q_index++){
                for (MKL_INT64 x_index = 0; x_index <= Nx_Gtilde - 1; x_index++){
                    for (MKL_INT64 alpha = 0; alpha <= Jj_dim - 1; alpha++){
                        for (MKL_INT64 alphaprime = 0; alphaprime <= Jj_dim - 1; alphaprime++){
                            Gtilde_store[alpha * Jj_dim * Np_3N * Nq_3N * Nx_Gtilde + alphaprime * Np_3N * Nq_3N * Nx_Gtilde + p_index * Nq_3N * Nx_Gtilde + q_index * Nx_Gtilde + x_index]
                                = Gtilde_new (p_3N[p_index], q_3N[q_index], x_Gtilde[x_index], alpha, alphaprime, Jj_dim, Lmax, L12_Jj, l3_Jj, Atilde_store, two_J);
                            
                            //counter += 1;
                            //frac_n = (100*counter)/fullsize;
                            //if (frac_n>frac_o){cout << frac_n << "%" << endl; frac_o=frac_n;}
                        }
                    }
                }
            }
        }
    }
    std::cout << "   - Done!" << std::endl;
    /* Free memory */
    delete [] Atilde_store;
    delete [] SixJ_array;
}

/* Indexing: S[j*Nx + i + k*Nx(Nx-1)] for S_i(x) where x is in bin j, and k is the coefficient we need
 * We have a function evaluation at x_i. 
 * This means we use x_j when x is in bin j, i.e. (x-x_j). */
double interpolate_using_spline_array(double* S_array, int Np, double* p_array, int idx_pi, double p, int idx_p){
    double interpolate = 0;

    double pi = p_array[idx_pi];

    if (p <= p_array[Np - 1]){
        interpolate = (
                                                           S_array[idx_pi * Np + idx_p]
                            + (p - pi)                   * S_array[idx_pi * Np + idx_p +     Np * (Np - 1)]
                            + (p - pi)*(p - pi)          * S_array[idx_pi * Np + idx_p + 2 * Np * (Np - 1)]
                            + (p - pi)*(p - pi)*(p - pi) * S_array[idx_pi * Np + idx_p + 3 * Np * (Np - 1)]
                        );
    }

    return interpolate;
}

void calculate_angular_quadrature_grids(double* x_array, double* wx_array, int Nx){
    calc_gauss_points (x_array, wx_array, -1., 1., Nx);
}

double extract_potential_element_from_array(int L, int Lp, int J, int S, bool coupled, double* V_array){
	if (J<abs(S-L) or J>S+L or J<abs(S-Lp) or J>S+Lp){
		std::cout << L <<" "<< Lp <<" "<< J <<" "<< S << std::endl;
		raise_error("Encountered unphysical state in LS-solver");
	}

    double potential_element = NAN;

    /* Check OPEP.cpp to see which elements in V_array correponds to which angular momenta */
    if (coupled){
        if (L==Lp and L<J){         // --
            potential_element = V_array[2];
        }
        else if (L==Lp and L>J){    // ++
            potential_element = V_array[5];
        }
        else if (L<Lp){             // -+
            potential_element = V_array[3];
        }
        else{                       // +-
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

/* Construct 2N potential matrices <k|v|k_p> for all 3N partial wave channels */
void calculate_potential_matrices_array(double* V_unco_array,
                                        double* V_coup_array,
                                        int Np, double* p_array, double* wp_array,
                                        int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array,
                                        potential_model* pot_ptr_nn,
                                        potential_model* pot_ptr_np){

    double V_IS_elements [6];	// Isoscalar (IS)
	double V_nn_elements [6];	// neutron-neutron (nn)
	double V_np_elements [6];	// neutron-proton (np)

    int Np1 = Np+1;

	double p_r=0, p_in=0, p_c=0, p_out=0;
    /* Row state */
    for (int idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
        int L_r = L_2N_array[idx_alpha_r];
        int S_r = S_2N_array[idx_alpha_r];
        int J_r = J_2N_array[idx_alpha_r];
        int T_r = T_2N_array[idx_alpha_r];

        /* Column state */
        for (int idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
            int L_c = L_2N_array[idx_alpha_c];
            int S_c = S_2N_array[idx_alpha_c];
            int J_c = J_2N_array[idx_alpha_c];
            int T_c = T_2N_array[idx_alpha_c];

            /* Check if possible channel through interaction */
            if (T_r==T_c and J_r==J_c and S_r==S_c and abs(L_r-L_c)<=2){

                /* Detemine if this is a coupled channel: NOTE THAT THE POTENTIAL MODEL
                 * REGISTERS 3P0 AS COUPLED ALTHOUGH IT'S AN UNCOUPLED MATRIX, HENCE coupled_model */
	            bool coupled_matrix = false;
                bool coupled_model  = false;
	            if (L_r!=L_c or (L_r==L_c and L_r!=J_r and J_r!=0)){
	            	coupled_matrix  = true;
	            }
                if (L_r!=L_c or (L_r==L_c and L_r!=J_r)){
	            	coupled_model   = true;
	            }

                /* Skip redundant calculations by only doing the coupled calculation when L_r<L_c */
                if (coupled_matrix){
                    if ( (L_r<L_c)==false ){
                        continue;
                    }
                }

                for (int idx_p_r=0; idx_p_r<Np; idx_p_r++){
                    p_r = p_array[idx_p_r];
                    p_in = p_r;//*hbarc;

                    for (int idx_p_c=0; idx_p_c<Np; idx_p_c++){
	                    p_c = p_array[idx_p_c];
                        p_out = p_c;//*hbarc;

	                    /* We create an isoscalar potential */
	                    if (T_r==1){ // Interaction can be either nn or np
                            pot_ptr_nn->V(p_in, p_out, coupled_model, S_r, J_r, T_r, V_nn_elements);
                            pot_ptr_np->V(p_in, p_out, coupled_model, S_r, J_r, T_r, V_np_elements);
	                    	for (int idx_element=0; idx_element<6; idx_element++){
	                    		V_IS_elements[idx_element] = (1./3)*V_np_elements[idx_element] + (2./3)*V_nn_elements[idx_element];
	                    	}
                        }
                        else{ 	   // Interaction must be np
                            pot_ptr_np->V(p_in, p_out, coupled_model, S_r, J_r, T_r, V_IS_elements);
                        }
                        //pot_ptr_np->V(p_in, p_out, coupled, S_r, J_r, T_r, V_IS_elements);

                        //double norm_fac = (M_PI/2);
                        double norm_fac = 1;
	                    
                        /* Write element to potential matrix V_array */
                        if (coupled_matrix){
                            int step_V_coup = J_r-1;
                            //std::cout << idx_p_r << " " << idx_p_r+Np1 << std::endl;
			            	V_coup_array[step_V_coup*4*Np1*Np1 +  idx_p_r       *2*Np1 + idx_p_c]	    = norm_fac*extract_potential_element_from_array(J_r-1, J_r-1, J_r, S_r, coupled_matrix, V_IS_elements);
			            	V_coup_array[step_V_coup*4*Np1*Np1 +  idx_p_r       *2*Np1 + idx_p_c + Np1] = -norm_fac*extract_potential_element_from_array(J_r-1, J_r+1, J_r, S_r, coupled_matrix, V_IS_elements);
			            	V_coup_array[step_V_coup*4*Np1*Np1 + (idx_p_r + Np1)*2*Np1 + idx_p_c]	    = -norm_fac*extract_potential_element_from_array(J_r+1, J_r-1, J_r, S_r, coupled_matrix, V_IS_elements);
			            	V_coup_array[step_V_coup*4*Np1*Np1 + (idx_p_r + Np1)*2*Np1 + idx_p_c + Np1] = norm_fac*extract_potential_element_from_array(J_r+1, J_r+1, J_r, S_r, coupled_matrix, V_IS_elements);
			            
                            //if (J_r==1){
                            //    std::cout << V_coup_array[step_V_coup*4*Np*Np +  idx_p_r      *2*Np + idx_p_c]	    << std::endl;
                            //    std::cout << V_coup_array[step_V_coup*4*Np*Np + (idx_p_r + Np)*2*Np + idx_p_c]	    << std::endl;
			            	//    std::cout << V_coup_array[step_V_coup*4*Np*Np +  idx_p_r      *2*Np + idx_p_c + Np] << std::endl;
			            	//    std::cout << V_coup_array[step_V_coup*4*Np*Np + (idx_p_r + Np)*2*Np + idx_p_c + Np] << std::endl;
                            //}
                        }
			            else{
                            int step_V_unco = 2*J_r + S_r;
			            	V_unco_array[step_V_unco*Np1*Np1 + idx_p_r*Np1 + idx_p_c] = norm_fac*extract_potential_element_from_array(L_r, L_c, J_r, S_r, coupled_matrix, V_IS_elements);

                            //if (S_r==0 and J_r==0){
                            //    std::cout << V_unco_array[step_V_unco*Np*Np + idx_p_r*Np + idx_p_c] << std::endl;
                            //}
                        }
                    }
                }
            }
        }
    }
}

void calculate_faddeev_convergence(double* state_array,
                                   double* P123_array, int Nalpha_P123, int Np_P123, int Nq_P123,
                                   int Np, double* p_array, double* wp_array,
                                   int Nq, double* q_array, double* wq_array,
                                   int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array, int* l_3N_array, int* two_j_3N_array,
                                   int two_T, int two_J, int PAR, int J_2N_max,
                                   potential_model* pot_ptr_nn,
                                   potential_model* pot_ptr_np){

    bool print_calculation_steps = true;

    /* pw basis size */
    int basis_size = Nalpha*Np*Nq;

    /* Potential terms arrays */
    int Np1 = Np+1;
    int V_unco_array_size = Np1*Np1   * 2*(J_2N_max+1);
    int V_coup_array_size = Np1*Np1*4 *   (J_2N_max+1);
    double* V_unco_array = new double [V_unco_array_size];
    double* V_coup_array = new double [V_coup_array_size];

    /* Ensure arrays contain zeros */
    for (int i=0; i<V_unco_array_size; i++){
        V_unco_array[i] = 0;
    }
    for (int i=0; i<V_coup_array_size; i++){
        V_coup_array[i] = 0;
    }

    /* Pre-calculating potentials */
    std::cout << "Pre-calculating potential matrices ... " << std::flush;
    calculate_potential_matrices_array(V_unco_array,
                                       V_coup_array,
                                       Np, p_array, wp_array,
                                       Nalpha, L_2N_array, S_2N_array, J_2N_array, T_2N_array,
                                       pot_ptr_nn,
                                       pot_ptr_np);
    std::cout << "Done." << std::endl;

    /* K-matrix */
    double* K_array = new double [basis_size * basis_size];
    
    /* Eigenvalues of K */
    double *lambda_array = new double [basis_size];

    /* Eigenvectors of K */
    double *psi_array = new double [basis_size * basis_size];

    /* Triton ground-state energy (to be determined) */
    double Z = -7, E=0, lambda=0, Z_step_length=0.5;

    // I found lambda=1.17705 for Z=-1.5382 for Jmax=2 (TOOK VERY LONG TIME)

    //if (Nq==10){
    //    if(J_2N_max==1){
    //        Z = -7.429;
    //    }
    //    else if (J_2N_max==2){
    //        Z = -7.592;
    //    }
    //}
    
    bool lambda_equals_one = false;
    bool Z_too_negative    = false;

    while (Z>=-9){
        if(print_calculation_steps){ std::cout << "Calculating K-matrix for Z=" << Z << " ... " << std::flush;}
        
        for (int i=0; i<basis_size*basis_size; i++){
            K_array[i] = 0;
        }

        iterate_faddeev(K_array, Z,
                        P123_array, Nalpha_P123, Np_P123, Nq_P123,
                        V_unco_array,
                        V_coup_array,
        	            Np, p_array, wp_array,
        	            Nq, q_array, wq_array,
        	            Nalpha, L_2N_array, S_2N_array, J_2N_array, T_2N_array, l_3N_array, two_j_3N_array);

        //construct_K_by_MM_multiplication(K_array, Z,
        //                                 P123_array, Nalpha_P123, Np_P123, Nq_P123,
        //                                 V_unco_array,
        //                                 V_coup_array,
        //	                             Np, p_array, wp_array,
        //	                             Nq, q_array, wq_array,
        //	                             Nalpha, L_2N_array, S_2N_array, J_2N_array, T_2N_array, l_3N_array, two_j_3N_array);            
                                
        if(print_calculation_steps){ std::cout << "Done." << std::endl;}
        
        /* Find eigenvalues lambda of K */
        double lambda_current = 0;
        if(print_calculation_steps){ std::cout << "Diagonalise K-matrix ... " << std::flush;}
        
        find_eigenvalues(K_array, lambda_array, psi_array, basis_size);
        for (int i=0; i<basis_size; i++){
            if (lambda_current < lambda_array[i]){
                lambda_current = lambda_array[i];
            }
        }
        
        if(print_calculation_steps){ std::cout << "Done." << std::endl;}
        
        /* See if we got closer to lambda=1 */
        if (lambda < lambda_current){
            lambda = lambda_current;
            E = Z;
        }

        printf("Highest lambda = %.5f for Z = %.4f. Highest lambda found so far is %.5f for E = %.4f. \n", lambda_current, Z, lambda, E);
        Z -= Z_step_length;
    }
}

void iterate_faddeev(double* K_array,
                     double Z,
                     double* P123_array, int Nalpha_P123, int Np_P123, int Nq_P123,
                     double* V_unco_array,
                     double* V_coup_array,
                     int Np, double* p_array, double* wp_array,
                     int Nq, double* q_array, double* wq_array,
                     int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array, int* l_3N_array, int* two_j_3N_array){

    int basis_size = Nalpha*Np*Nq;

    /* P123-matrix element of interest */
    double P123_element = 0;
    /* T-matrix element of interest */
    double t_element = 0;

    /* t-matrix arrays. These are fully rewritten on a call to calculate_t_element */
    cfloatType* t_unco_array = new cfloatType [  (Np+1)*(Np+1)];
    cfloatType* t_coup_array = new cfloatType [4*(Np+1)*(Np+1)];

    /* K-matrix summation terms for each column-state */
    double* K_summation_terms = new double [basis_size];
    for (int idx=0; idx<basis_size; idx++){
        K_summation_terms[idx] = 0;
    }

    //int exl_array [6] = {0,1,2,3,6,7};
    //int exl_idx1 = 2;
    //int exl_idx2 = 3;
    //int exl_idx3 = 6;
    //int exl_idx4 = 7;
    //std::cout << "\n Ignoring indices: " << exl_idx1 << " and " << exl_idx2 << std::endl;
    /* Row state */
	for (int idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){

        //if (idx_alpha_r==exl_idx1 or idx_alpha_r==exl_idx2 or idx_alpha_r==exl_idx3 or idx_alpha_r==exl_idx4){continue;}

        /* Projection state quantum numbers (rows) */
        int L_r     = L_2N_array[idx_alpha_r];
        int S_r     = S_2N_array[idx_alpha_r];
        int J_r     = J_2N_array[idx_alpha_r];
        int T_r     = T_2N_array[idx_alpha_r];
        int two_l_r = l_3N_array[idx_alpha_r];
        int two_j_r = two_j_3N_array[idx_alpha_r];
        
        /* Loop over alpha prime summation */
        for (int idx_alpha_p=0; idx_alpha_p<Nalpha; idx_alpha_p++){

            //if (idx_alpha_p==exl_idx1 or idx_alpha_p==exl_idx2 or idx_alpha_p==exl_idx3 or idx_alpha_p==exl_idx4){continue;}

            int L_p     = L_2N_array[idx_alpha_p];
            int S_p     = S_2N_array[idx_alpha_p];
            int J_p     = J_2N_array[idx_alpha_p];
            int T_p     = T_2N_array[idx_alpha_p];
            int two_l_p = l_3N_array[idx_alpha_p];
            int two_j_p = two_j_3N_array[idx_alpha_p];

            /* The two-body force can only allow for L!=L_pp, everything else must be the same */
            if (two_l_r==two_l_p and two_j_r==two_j_p and T_r==T_p and J_r==J_p and S_r==S_p and abs(L_r-L_p)<=2){
                bool coupled = L_r!=L_p or (L_r==L_p and L_r!=J_p and J_r!=0);

                double*     V_mat_array_ptr = NULL;
                cfloatType* t_mat_array_ptr = NULL;
	            if (coupled){ // coupled interaction
                    int step_V_coup = J_r-1;
	            	V_mat_array_ptr = &V_coup_array[step_V_coup*4*(Np+1)*(Np+1)];
                    t_mat_array_ptr = t_coup_array;
	            }
                else{ // uncoupled
                    int step_V_unco = 2*J_r + S_r;
                    V_mat_array_ptr = &V_unco_array[step_V_unco*(Np+1)*(Np+1)];
                    t_mat_array_ptr = t_unco_array;
                }

                for (int idx_q_r=0; idx_q_r<Nq; idx_q_r++){
                    double q_kin_term = 0.75*q_array[idx_q_r]*q_array[idx_q_r]/MN;

                    /* Energy for LS-solver */
                    double E_LS = Z - q_kin_term;
                    /* Solve Lippmann-Schwinger equation */
                    calculate_t_element(V_mat_array_ptr,
                                        t_mat_array_ptr,
                                        coupled,
		    	    		            E_LS, MN,
		    	    		            Np, p_array, wp_array);

                    for (int idx_p_r=0; idx_p_r<Np; idx_p_r++){
                        double p_kin_term = p_array[idx_p_r]*p_array[idx_p_r]/MN;

                        int idx_K_row = idx_alpha_r*Nq*Np + idx_q_r*Np + idx_p_r;
                                                        
                        /* Loop over integration momenta */
                        for (int idx_p_p=0; idx_p_p<Np; idx_p_p++){
                            /* Extract correct matrix element */
                            if (coupled){
                                if(L_r<L_p){        // Upper right
                                    t_element = t_mat_array_ptr[ idx_p_r      *2*(Np+1) + idx_p_p + Np+1].real();
                                }
                                else if (L_r>L_p){  // Lower left
                                    t_element = t_mat_array_ptr[(idx_p_r+Np+1)*2*(Np+1) + idx_p_p].real();
                                }
                                else if (L_r<J_r){  // Upper left
                                    t_element = t_mat_array_ptr[ idx_p_r      *2*(Np+1) + idx_p_p].real();
                                }
                                else{               // Lower right
                                    t_element = t_mat_array_ptr[(idx_p_r+Np+1)*2*(Np+1) + idx_p_p + Np+1].real();
                                }
                            }
                            else{
                                t_element = t_mat_array_ptr[idx_p_r*(Np+1) + idx_p_p].real();
                            }

                            /* Column state */
                            for (int idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
                                for (int idx_q_c=0; idx_q_c<Nq; idx_q_c++){
                                    for (int idx_p_c=0; idx_p_c<Np; idx_p_c++){
                                        int idx_K_col = idx_alpha_c*Nq*Np + idx_q_c*Np + idx_p_c;

                                        /* Retrieve P123 from pre-calculated P123_array - this INCLUDES spline-functionality*/ 
                                        //int P123_idx = (idx_alpha_c*Nq_P123*Np_P123 + idx_q_c*Np_P123 + idx_p_c)*Np_P123*Nq_P123*Nalpha_P123
                                        //              + idx_alpha_p*Nq_P123*Np_P123 + idx_q_r*Np_P123 + idx_p_p;
                                        int P123_idx = (idx_alpha_p*Nq_P123*Np_P123 + idx_q_r*Np_P123 + idx_p_p)*Np_P123*Nq_P123*Nalpha_P123
                                                      + idx_alpha_c*Nq_P123*Np_P123 + idx_q_c*Np_P123 + idx_p_c;

                                        P123_element = P123_array[P123_idx];

                                        K_summation_terms[idx_K_col] += wp_array[idx_p_p] * p_array[idx_p_p] * p_array[idx_p_p] * t_element * P123_element;
                                    }
                                }
                            }
                        }

                        /* 3-nucleon free Green's function */
                        double G0 = 1./(Z - p_kin_term - q_kin_term );
        
                        for (int idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
                            
                            //if (idx_alpha_c==exl_idx1 or idx_alpha_c==exl_idx2 or idx_alpha_c==exl_idx3 or idx_alpha_c==exl_idx4){continue;}

                            for (int idx_q_c=0; idx_q_c<Nq; idx_q_c++){
                                for (int idx_p_c=0; idx_p_c<Np; idx_p_c++){
                                    int idx_K_col = idx_alpha_c*Nq*Np + idx_q_c*Np + idx_p_c;
                                    int idx_K = idx_K_row*basis_size + idx_K_col;
                                    K_array[idx_K] += 2*G0*K_summation_terms[idx_K_col];
                                }
                            }
                        }
        
                        /* Reset K-summation terms */
                        for (int idx=0; idx<basis_size; idx++){
                            K_summation_terms[idx] = 0;
                        }
                    }
                }
            }
        }
    }

    delete [] t_unco_array;
    delete [] t_coup_array;
}


void dot_MM(float *A, float *B, float *C, int N, int K, int M){
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, M, K, 1.0, A, K, B, M, 0.0, C, M);
}
void dot_MM(double *A, double *B, double *C, int N, int K, int M){
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, M, K, 1.0, A, K, B, M, 0.0, C, M);
}
void construct_K_by_MM_multiplication(double* K_array,
                                      double Z,
                                      double* P123_array, int Nalpha_P123, int Np_P123, int Nq_P123,
                                      double* V_unco_array,
                                      double* V_coup_array,
                                      int Np, double* p_array, double* wp_array,
                                      int Nq, double* q_array, double* wq_array,
                                      int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array, int* l_3N_array, int* two_j_3N_array){
    
    int K_dim    = Nalpha*Np*Nq;
    int K_dim_sq = K_dim*K_dim;
    double* G0_array = new double [K_dim_sq];
    double* T_array  = new double [K_dim_sq];

    /* Set all elements to zero */
    for (int i=0; i<K_dim_sq; i++){
        G0_array[i] = 0;
        T_array[i]  = 0;
    }

    int idx = 0;
    double* V_mat_array_ptr;
    cfloatType* t_unco_array = new cfloatType [  (Np+1)*(Np+1)];
    cfloatType* t_coup_array = new cfloatType [4*(Np+1)*(Np+1)];
    /* Bra and ket summation for the three matrices */
    for (int idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
        for (int idx_q_r=0; idx_q_r<Nq; idx_q_r++){
            for (int idx_p_r=0; idx_p_r<Np; idx_p_r++){
                
                int L_r     = L_2N_array[idx_alpha_r];
                int S_r     = S_2N_array[idx_alpha_r];
                int J_r     = J_2N_array[idx_alpha_r];
                int T_r     = T_2N_array[idx_alpha_r];
                int two_l_r = l_3N_array[idx_alpha_r];
                int two_j_r = two_j_3N_array[idx_alpha_r];

                double p = p_array[idx_p_r];
                double q = q_array[idx_q_r];
                double p_kinetic_term =      p*p/MN;
                double q_kinetic_term = 0.75*q*q/MN;

                /* G0-matrix construction */
                idx = (idx_alpha_r*Nq*Np + idx_q_r*Np + idx_p_r)*Np*Nq*Nalpha
                     + idx_alpha_r*Nq*Np + idx_q_r*Np + idx_p_r;
                G0_array[idx] = 1./(Z - p_kinetic_term - q_kinetic_term);

                for (int idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
                    for (int idx_q_c=0; idx_q_c<Nq; idx_q_c++){
                        for (int idx_p_c=0; idx_p_c<Np; idx_p_c++){

                            int L_c     = L_2N_array[idx_alpha_c];
                            int S_c     = S_2N_array[idx_alpha_c];
                            int J_c     = J_2N_array[idx_alpha_c];
                            int T_c     = T_2N_array[idx_alpha_c];
                            int two_l_c = l_3N_array[idx_alpha_c];
                            int two_j_c = two_j_3N_array[idx_alpha_c];

                            if (idx_q_c==idx_q_r and T_r==T_c and S_r==S_c and J_r==J_c and two_l_r==two_l_c and two_j_r==two_j_c and abs(L_r-L_c)<=2){

                                bool coupled = (L_r!=L_c or (L_r==L_c and L_r!=J_r and J_r!=0));
                                cfloatType* t_mat_array_ptr = NULL;
                                if (coupled){
                                    V_mat_array_ptr = &V_coup_array[(J_r-1)*4*(Np+1)*(Np+1)];
                                    t_mat_array_ptr = t_coup_array;
                                    for (int i=0; i<4*(Np+1)*(Np+1); i++){
                                        t_mat_array_ptr[i] = 0;
                                    }
                                }
                                else{
                                    V_mat_array_ptr = &V_unco_array[(2*J_r+S_r)*(Np+1)*(Np+1)];
                                    t_mat_array_ptr = t_unco_array;
                                    for (int i=0; i<(Np+1)*(Np+1); i++){
                                        t_mat_array_ptr[i] = 0;
                                    }
                                }

                                double E_LS = Z - q_kinetic_term;
                                calculate_t_element(V_mat_array_ptr,
                                                    t_mat_array_ptr,
                                                    coupled,
		    	    		                        E_LS, MN,
		    	    		                        Np, p_array, wp_array);
                                
                                double t_element = 0;
                                if (coupled){
                                    if(L_r<L_c){        // Upper right
                                        t_element = t_mat_array_ptr[ idx_p_r      *2*(Np+1) + idx_p_c + Np+1].real();
                                    }
                                    else if (L_r>L_c){  // Lower left
                                        t_element = t_mat_array_ptr[(idx_p_r+Np+1)*2*(Np+1) + idx_p_c].real();
                                    }
                                    else if (L_r<J_r){  // Upper left
                                        t_element = t_mat_array_ptr[ idx_p_r      *2*(Np+1) + idx_p_c].real();
                                    }
                                    else{               // Lower right
                                        t_element = t_mat_array_ptr[(idx_p_r+Np+1)*2*(Np+1) + idx_p_c + Np+1].real();
                                    }
                                }
                                else{
                                    t_element = t_mat_array_ptr[idx_p_r*(Np+1) + idx_p_c].real();
                                }

                                idx = (idx_alpha_r*Nq*Np + idx_q_r*Np + idx_p_r)*Np*Nq*Nalpha
                                     + idx_alpha_c*Nq*Np + idx_q_c*Np + idx_p_c;
                                
                                T_array[idx] = wp_array[idx_p_c]*p_array[idx_p_c]*p_array[idx_p_c]*t_element;
                            }
                        }
                    }
                }
            }
        }
    }

    double* temp_array = new double [K_dim_sq];

    for (int i=0; i<K_dim_sq; i++){
        temp_array[i] = 0;
    }

    dot_MM(G0_array, T_array, temp_array, K_dim, K_dim, K_dim);

    dot_MM(temp_array, P123_array, K_array, K_dim, K_dim, K_dim);

    delete [] t_unco_array;
    delete [] t_coup_array;
    delete [] G0_array;
    delete [] T_array;
    delete [] temp_array;
}
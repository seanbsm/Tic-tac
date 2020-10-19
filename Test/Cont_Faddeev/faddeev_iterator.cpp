
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
    int ldvl = 0;   // leading dimension of vl

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

void calculate_faddeev_convergence(double* state_array,
                                   double* P123_array,
                                   int Np, double* p_array, double* wp_array,
                                   int Nq, double* q_array, double* wq_array,
                                   int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array, int* l_3N_array, int* two_j_3N_array,
                                   int two_T, int two_J, int PAR,
                                   potential_model* pot_ptr_nn,
                                   potential_model* pot_ptr_np){
    
    /* pw basis size */
    int basis_size = Nalpha*Np*Nq;

    /* K-matrix */
    double* K_array = new double [basis_size * basis_size];
    
    /* Eigenvalues of K */
    double *lambda_array = new double [basis_size];

    /* Eigenvectors of K */
    double *psi_array = new double [basis_size * basis_size];

    /* K(z) value) */
    double K = 0;

    /* Triton ground-state energy (to be determined) */
    double Z = 0;

    cout << "Calculating K-matrix" << endl;
    /* Row state */
	for (int idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
        for (int idx_q_r=0; idx_q_r<Nq; idx_q_r++){
            for (int idx_p_r=0; idx_p_r<Np; idx_p_r++){
                int idx_K_row = idx_alpha_r*Nq*Np + idx_q_r*Np + idx_p_r;

                /* Column state */
                for (int idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
                    for (int idx_q_c=0; idx_q_c<Nq; idx_q_c++){
                        for (int idx_p_c=0; idx_p_c<Np; idx_p_c++){
                            int idx_K_col = idx_alpha_c*Nq*Np + idx_q_c*Np + idx_p_c;

                            K = iterate_faddeev(Z,
                                                P123_array,
    		                	                Np, p_array, wp_array,
    		                	                Nq, q_array, wq_array,
    		                	                Nalpha, L_2N_array, S_2N_array, J_2N_array, T_2N_array, l_3N_array, two_j_3N_array,
    		                	                idx_alpha_r, idx_p_r, idx_q_r,
                                                idx_alpha_c, idx_p_c, idx_q_c,
    		                	                two_T, two_J, PAR,
    		                	                pot_ptr_nn,
    		                	                pot_ptr_np);

                            int idx_K = idx_K_row*basis_size + idx_K_col;
                            K_array[idx_K] = K;
                        }
                    }
                }
            }
        }
    }

    find_eigenvalues(K_array, lambda_array, psi_array, basis_size);

}

double iterate_faddeev(double Z,
                       double* P123_array,
                       int Np, double* p_array, double* wp_array,
                       int Nq, double* q_array, double* wq_array,
                       int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array, int* l_3N_array, int* two_j_3N_array,
                       int idx_alpha_r, int idx_p_r, int idx_q_r,
                       int idx_alpha_c, int idx_p_c, int idx_q_c,
                       int two_T, int two_J, int PAR,
                       potential_model* pot_ptr_nn,
                       potential_model* pot_ptr_np){

    int idx_alpha = idx_alpha_r;
    int idx_p = idx_p_r;
    int idx_q = idx_q_r;

    if (idx_q_r != idx_q_c){
        return 0;
    }

    /* Projection state quantum numbers (alpha) */
    int L     = L_2N_array[idx_alpha];
    int S     = S_2N_array[idx_alpha];
    int J     = J_2N_array[idx_alpha];
    int T     = T_2N_array[idx_alpha];
    int two_l = l_3N_array[idx_alpha];
    int two_j = two_j_3N_array[idx_alpha];

    /* Upper-case is for the 2N-pair, lower-case is for orbiting spectator */
    int L_p  = 0, S_p  = 0, J_p  = 0, T_p  = 0, two_l_p  = 0, two_j_p  = 0;
    int L_pp = 0, S_pp = 0, J_pp = 0, T_pp = 0, two_l_pp = 0, two_j_pp = 0;
    /* P123-matrix element of interest */
    double P123_element = 0;
    /* T-matrix element of interest */
    double t_element = 0;
    /* Faddeev summation term */
    double total_sum = 0;
    /* Loop over alpha prime summation */
    for (int idx_alpha_p=0; idx_alpha_p<Nalpha; idx_alpha_p++){
        std::cout << idx_alpha_p << " / " << Nalpha << std::endl;
        L_p     = L_2N_array[idx_alpha_p];
        S_p     = S_2N_array[idx_alpha_p];
        J_p     = J_2N_array[idx_alpha_p];
        T_p     = T_2N_array[idx_alpha_p];
        two_l_p = l_3N_array[idx_alpha_p];
        two_j_p = two_j_3N_array[idx_alpha_p];
        
        /* Loop over alpha double-prime summation */
        for (int idx_alpha_pp=0; idx_alpha_pp<Nalpha; idx_alpha_pp++){
            L_pp     = L_2N_array[idx_alpha_pp];
            S_pp     = S_2N_array[idx_alpha_pp];
            J_pp     = J_2N_array[idx_alpha_pp];
            T_pp     = T_2N_array[idx_alpha_pp];
            two_l_pp = l_3N_array[idx_alpha_pp];
            two_j_pp = two_j_3N_array[idx_alpha_pp];

            /* The two-body force can only allow for L!=L_pp, everything else must be the same */
            if (two_l==two_l_pp and two_j==two_j_pp and T==T_pp and J==J_pp and S==S_pp and abs(L-L_pp)<=2){
                if (idx_alpha_pp == idx_alpha_c){
                    /* Loop over p prime integral */
                    for (int idx_p_p=0; idx_p_p<Np; idx_p_p++){
                        if (idx_p_p == idx_p_c){

                            /* Solve Lippmann-Schwinger equation for current q', q, alpha', and alpha */
                            /* !!! NOTE I USE NUCLEON MASS HERE, THIS IS SLIGHTLY WRONG BUT SHOULD MATTER LITTLE !!! */
                            double E_LS = Z - q_array[idx_q]*q_array[idx_q] / (2 * 3*MN/4);
                            double q_LS = sqrt(E_LS*MN);
                            t_element = calculate_t_element(L, L_pp, S, J, T,
			                			                    E_LS, q_LS, MN,
			                			                    Np, p_array, wp_array,
			                			                    idx_p, idx_p_p,
			                			                    pot_ptr_nn, pot_ptr_np);
    
                            /* Retrieve P123 from pre-calculated P123_array - this INCLUDES spline-functionality*/
                            int P123_idx = (idx_alpha*Nq*Np + idx_q*Np + idx_p)*Np*Nq*Nalpha + idx_alpha_pp*Nq*Np + idx_q*Np + idx_p_p;
                            P123_element = P123_array[P123_idx];
    
                            /* Evaluate integral term and add to total sum */
                            total_sum +=  wp_array[idx_p_p] * p_array[idx_p_p] * p_array[idx_p_p]
                                        * t_element * P123_element;
                        }
                    }
                }
            }
        }
    }

    /* First faddeev state summation end */

    /* 3-nucleon free Green's function */
    double mu1 = 0.5*MN;
    double mu2 = 2.0*MN/3.;
    double G0 = 1./(Z - 0.5*p_array[idx_p]*p_array[idx_p]/mu1 - 0.5*q_array[idx_q]*q_array[idx_q]/mu2 );

    return total_sum*G0;
}
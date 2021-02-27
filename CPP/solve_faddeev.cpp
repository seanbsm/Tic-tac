
#include "solve_faddeev.h"

void calculate_CPVC_product_row(){
}

void pade_method_solve(){
}

void direct_sparse_solve(){
}

/* Solves the Faddeev equations
 * U = P*V + P*V*G*U
 * on the form L*U = R, where L and R are the left-
 * and right-handed sides of the equations, given by 
 * L = 1 - P*V*G
 * R = P*V
 * Since G is expressed in an SWP basis, we also must include the basis-transormation matrices C */
void direct_dense_solve(cdouble* U_array,
                        cdouble* G_array,
                        double* P123_array,
                        double* C_WP_unco_array,
				        double* C_WP_coup_array,
				        double* V_WP_unco_array,
                        double* V_WP_coup_array,
                        int Nq_WP,
				        int Np_WP,
				        int Nalpha,
                        int* L_2N_array,
                        int* S_2N_array,
                        int* J_2N_array,
                        int* T_2N_array,
                        int* L_1N_array, 
                        int* two_J_1N_array,
                        int* two_J_3N_array,
                        int* two_T_3N_array){
    
    int mat_dim    = Nalpha * Nq_WP * Np_WP;
    int mat_dim_sq = mat_dim * mat_dim;

    cdouble* L_array = new cdouble [mat_dim_sq];

    /* Temporary loop variables */
    int mat_2N_dim = 0;
    int chn_2N_idx = 0;
    double* V_array_ptr = NULL;
    double* C_array_ptr = NULL;

    auto setup_start = std::chrono::system_clock::now();
    /* <X_i'j'^alpha'| - loops */
    for (int idx_alpha_r = 0; idx_alpha_r < Nalpha; idx_alpha_r++){
        for (int idx_p_bin_r = 0; idx_p_bin_r < Np_WP; idx_p_bin_r++){
            for (int idx_q_bin_r = 0; idx_q_bin_r < Nq_WP; idx_q_bin_r++){

                /* |X_ij^alpha> - loops */
                for (int idx_alpha_c = 0; idx_alpha_c < Nalpha; idx_alpha_c++){

                    int L_2N_c     = L_2N_array[idx_alpha_c];
                    int S_2N_c     = S_2N_array[idx_alpha_c];
                    int J_2N_c     = J_2N_array[idx_alpha_c];
                    int T_2N_c     = T_2N_array[idx_alpha_c];
                    int L_1N_c     = L_1N_array[idx_alpha_c];
                    int two_J_1N_c = two_J_1N_array[idx_alpha_c];
                    int two_J_3N_c = two_J_3N_array[idx_alpha_c];
                    int two_T_3N_c = two_T_3N_array[idx_alpha_c];

                    for (int idx_p_bin_c = 0; idx_p_bin_c < Np_WP; idx_p_bin_c++){
                        for (int idx_q_bin_c = 0; idx_q_bin_c < Nq_WP; idx_q_bin_c++){
                            
                            /* Index of G diagonal - equivalent to indexing of |X_ij^alpha> */
                            int G_idx = idx_alpha_c*Nq_WP*Np_WP + idx_q_bin_c*Np_WP + idx_p_bin_c;

                            double inner_product_PV    =  0;
                            cdouble inner_product_CGCT = {0, 0};

                            for (int idx_alpha_s = 0; idx_alpha_s < Nalpha; idx_alpha_s++){
                                
                                int L_2N_s     = L_2N_array[idx_alpha_s];
                                int S_2N_s     = S_2N_array[idx_alpha_s];
                                int J_2N_s     = J_2N_array[idx_alpha_s];
                                int T_2N_s     = T_2N_array[idx_alpha_s];
                                int L_1N_s     = L_1N_array[idx_alpha_s];
                                int two_J_1N_s = two_J_1N_array[idx_alpha_s];
                                int two_J_3N_s = two_J_3N_array[idx_alpha_s];
                                int two_T_3N_s = two_T_3N_array[idx_alpha_s];

                                int L_s  = L_2N_array[idx_alpha_s];
                                int S_s  = S_2N_array[idx_alpha_s];
                                int J_s  = J_2N_array[idx_alpha_s];
                                int T_s  = T_2N_array[idx_alpha_s];
                                int l_s  = L_1N_array[idx_alpha_s];
                                int j2_s = two_J_1N_array[idx_alpha_s];

                                /* Check if possible channel through interaction (only L_2N can vary through tensor-force) */
                                if ( T_2N_s==T_2N_c         &&
                                     J_2N_s==J_2N_c         &&
                                     S_2N_s==S_2N_c         &&
                                     abs(L_2N_s-L_2N_c)<=2  &&
                                     L_1N_s==L_1N_c         &&
                                     two_J_1N_s==two_J_1N_c &&
                                     two_J_3N_s==two_J_3N_c &&
                                     two_T_3N_s==two_T_3N_c ){

                                    /* Detemine if this is a coupled channel */
	                                if (L_2N_c!=J_2N_c and J_2N_c!=0){ // This counts 3P0 as uncoupled (which is intended)
	                                	mat_2N_dim = 2*Np_WP;
                                        chn_2N_idx = J_2N_s-1;
                                        V_array_ptr = &V_WP_coup_array[chn_2N_idx * mat_2N_dim * mat_2N_dim];
                                        C_array_ptr = &C_WP_coup_array[chn_2N_idx * mat_2N_dim * mat_2N_dim];
	                                }
                                    else{
                                        mat_2N_dim = Np_WP;
                                        chn_2N_idx = 2*J_2N_s + S_2N_s;
                                        V_array_ptr = &V_WP_unco_array[chn_2N_idx * mat_2N_dim * mat_2N_dim];
                                        C_array_ptr = &C_WP_unco_array[chn_2N_idx * mat_2N_dim * mat_2N_dim];
                                    }

                                    for (int idx_p_bin_s = 0; idx_p_bin_s < Np_WP; idx_p_bin_s++){
                                        for (int idx_q_bin_s = 0; idx_q_bin_s < Nq_WP; idx_q_bin_s++){

                                            /* Inner product P*V, to be appended to R_array */
                                            int V_idx = idx_p_bin_s * mat_2N_dim + idx_p_bin_c;
                                            int P_idx = (  idx_alpha_r*Nq_WP*Np_WP + idx_q_bin_r*Np_WP + idx_p_bin_r ) * mat_dim
                                                         + idx_alpha_s*Nq_WP*Np_WP + idx_q_bin_s*Np_WP + idx_p_bin_s;
                                            inner_product_PV += P123_array[P_idx] * V_array_ptr[V_idx];

                                            /* Inner product C*G*C^T, to be appended to L_array */
                                            int C_idx  = idx_p_bin_r * mat_2N_dim  + idx_p_bin_s;
                                            int CT_idx = idx_p_bin_c * mat_2N_dim  + idx_p_bin_s;
                                            int G_idx  = idx_alpha_s * Nq_WP*Np_WP + idx_q_bin_s*Np_WP + idx_p_bin_s;
                                            inner_product_CGCT += C_array_ptr[C_idx] * G_array[G_idx] * C_array_ptr[CT_idx];
                                        }
                                    }
                                }
                            }

                            /* Index of R and L */
                            int LR_idx = (  idx_alpha_c*Nq_WP*Np_WP + idx_q_bin_c*Np_WP + idx_p_bin_c ) * mat_dim
                                          + idx_alpha_r*Nq_WP*Np_WP + idx_q_bin_r*Np_WP + idx_p_bin_r;
                            
                            L_array[LR_idx] = inner_product_CGCT;
                            U_array[LR_idx] = {inner_product_PV, 0};
                        }
                    }
                }
            }
        }
    }
    auto setup_end = std::chrono::system_clock::now();
    std::chrono::duration<double> time_used = setup_end - setup_start;
	printf(" - Time used on setting up Faddeev equation:   %*.6f s \n", 10,time_used.count());

    /* Solve Faddeev (find U-array using L and R) */
    auto solve_start = std::chrono::system_clock::now();
    solve_MM(L_array, U_array, mat_dim);
    auto solve_end = std::chrono::system_clock::now();
    time_used = solve_end - solve_start;
    printf(" - Time used on solving Faddeev equation:      %*.6f s \n", 10,time_used.count());
}

void calculate_PVC_col(double*  col_array,
                       int      idx_alpha_c, int idx_p_c, int idx_q_c,
                       int      Nalpha,      int Nq_WP,   int Np_WP,
                       double** VC_CM_array,
                       double*  P123_val_array,
                       int*     P123_row_array,
                       int*     P123_col_array,
                       int      P123_dim){

    double* VC_subarray     = NULL;
    double* VC_subarray_col = NULL;

    int dense_dim   = Nalpha*Np_WP*Nq_WP;

    for (int idx_alpha_j=0; idx_alpha_j<Nalpha; idx_alpha_j++){
        
        int idx_VC_2N_block = idx_alpha_j*Nalpha + idx_alpha_c;
        VC_subarray = VC_CM_array[idx_VC_2N_block];

        /* Only do inner-product if VC is not zero due to conservation laws */
        if (VC_subarray!=NULL){

            VC_subarray_col = &VC_subarray[idx_p_c*Np_WP];
            
            /* Loop over rows of col-array */
            for (int idx_i=0; idx_i<dense_dim; idx_i++){
                
                /* CSR-format indexing */
                int idx_j_lower = P123_row_array[idx_i    ];
                int idx_j_upper = P123_row_array[idx_i + 1];
                
                /* Loop over inner-product indices (columns of P123)  */
                double inner_product_PVC = 0;
                for (int idx_j=idx_j_lower; idx_j<idx_j_upper; idx_j++){
                    int col_idx = P123_col_array[idx_j];
                    
                    /* Retrieve alpha_j index by using that
                     * the step-length per idx_alpha is Np_WP*Nq_WP */
                    div_t divresult1 = std::div(col_idx, Np_WP*Nq_WP);

                    /* Retrieve q_j and p_j index by using that
                     * the step-length per idx_q is Np_WP */
                    div_t divresult2 = std::div(divresult1.rem, Np_WP);

                    int idx_alpha_j_check = divresult1.quot;
                    int idx_q_j_check     = divresult2.quot;
                    if (idx_alpha_j==idx_alpha_j_check and idx_q_c==idx_q_j_check){

                        /* Extract p-momentum index from remainder of idx_q_j/Np_WP */
                        int idx_p_j = divresult2.rem;

                        /* Access arrays, this whole function is written to minimize these two calls */
                        double P_element  = P123_val_array[idx_j];
                        double VC_element = 0;//VC_subarray_col[idx_p_j];
                        inner_product_PVC += P_element * VC_element;
                    }
                }
                
                /* Write inner product to col_array of PVC-product */
                col_array[idx_i] = inner_product_PVC;
            }
        }
    }
}

void calculate_CPVC_col(double*  col_array,
                        int      idx_alpha_c, int idx_p_c, int idx_q_c,
                        int      Nalpha,      int Nq_WP,   int Np_WP,
                        double** CT_RM_array,
                        double** VC_CM_array,
                        double*  P123_val_array,
                        int*     P123_row_array,
                        int*     P123_col_array,
                        int      P123_dim){
    
    /* Generate PVC-column */
    double PVC_col [Nalpha*Nq_WP*Np_WP];
    /* Ensure PVC_col contains only zeroes */
    for (int idx=0; idx<Nalpha*Nq_WP*Np_WP; idx++){
        PVC_col[idx] = 0;
    }
    
    calculate_PVC_col(PVC_col,
                      idx_alpha_c, idx_p_c, idx_q_c,
                      Nalpha,      Nq_WP,   Np_WP,
                      VC_CM_array,
                      P123_val_array,
                      P123_row_array,
                      P123_col_array,
                      P123_dim);

    /* Generate (C^T x PVC)-column */
    double* CT_subarray     = NULL;
    double* CT_subarray_row = NULL;

    /* Loop over rows of col_array */
    for (int idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
        for (int idx_q_r=0; idx_q_r<Nq_WP; idx_q_r++){
            for (int idx_p_r=0; idx_p_r<Np_WP; idx_p_r++){

                double inner_product_CPVC = 0;

                /* Beginning of inner-product loops (index "i") */
                for (int idx_alpha_i=0; idx_alpha_i<Nalpha; idx_alpha_i++){
                    int idx_CT_2N_block = idx_alpha_r*Nalpha + idx_alpha_i;
                    CT_subarray = CT_RM_array[idx_CT_2N_block];

                    /* Only do inner-product if CT is not zero due to conservation laws */
                    if (CT_subarray!=NULL){
                        CT_subarray_row = &CT_subarray[idx_p_r*Np_WP];
                    
                        for (int idx_p_i=0; idx_p_i<Np_WP; idx_p_i++){
                            
                            int    idx_PVC     = idx_alpha_i*Nq_WP*Np_WP + idx_q_r*Np_WP + idx_p_i;
                            double PVC_element = P123_val_array[idx_PVC];

                            /* I'm not sure if this is the fastest ordering of the loops */
                            double CT_element  = 0;//CT_subarray_row[idx_p_i];

                            inner_product_CPVC += CT_element * PVC_element;
                        }
                    }
                }

                int idx_CPVC = idx_alpha_r*Nq_WP*Np_WP + idx_q_r*Np_WP + idx_p_r;

                col_array[idx_CPVC] = inner_product_CPVC;
            }
        }
    }
}

void solve_faddeev_equations(cdouble*  U_array,
                             cdouble** G_array,
                             double**  P123_sparse_ptr_val_array,
                             int**     P123_sparse_ptr_row_array,
                             int**     P123_sparse_ptr_col_array,
                             int*      P123_sparse_ptr_dim_array,
                             double*   C_WP_unco_array,
					         double*   C_WP_coup_array,
					         double*   V_WP_unco_array,
                             double*   V_WP_coup_array,
							 int       N_chn_3N,
    	                     int*      chn_3N_idx_array,
                             int       J_2N_max,
                             int       Nq_WP,
					         int       Np_WP,
					         int       Nalpha,
                             int*      L_2N_array,
                             int*      S_2N_array,
                             int*      J_2N_array,
                             int*      T_2N_array,
                             int*      L_1N_array, 
                             int*      two_J_1N_array,
                             int*      two_J_3N_array,
                             int*      two_T_3N_array){
    
    /* Create C^T-product pointer-arrays in row-major format */
    double*** CT_RM_ptr_array = new double** [N_chn_3N];
    create_CT_row_maj_3N_pointer_array(CT_RM_ptr_array,
                                       C_WP_unco_array,
					                   C_WP_coup_array,
                                       Np_WP,
                                       J_2N_max,
							           N_chn_3N,
    	                               chn_3N_idx_array,
					                   Nalpha,
                                       L_2N_array,
                                       S_2N_array,
                                       J_2N_array,
                                       T_2N_array);

    /* Create VC-product pointer-arrays in column-major format */
    double*** VC_CM_ptr_array  = new double** [N_chn_3N];
    create_VC_col_maj_3N_pointer_array(VC_CM_ptr_array,
                                       C_WP_unco_array,
					                   C_WP_coup_array,
					                   V_WP_unco_array,
                                       V_WP_coup_array,
					                   Np_WP,
                                       J_2N_max,
							           N_chn_3N,
    	                               chn_3N_idx_array,
					                   Nalpha,
                                       L_2N_array,
                                       S_2N_array,
                                       J_2N_array,
                                       T_2N_array);

    /* 3N-channel loop */
    for (int chn_3N=0; chn_3N< N_chn_3N; chn_3N++){

        int idx_alpha_lower = chn_3N_idx_array[chn_3N];
        int idx_alpha_upper = chn_3N_idx_array[chn_3N+1];
        int Nalpha_block    = idx_alpha_upper - idx_alpha_lower;

        int block_dense_dim = Nalpha_block * Nq_WP * Np_WP;

        /* P123 channel blocks */
        double* P123_val_array = P123_sparse_ptr_val_array[chn_3N];
        int*    P123_row_array = P123_sparse_ptr_row_array[chn_3N];
        int*    P123_col_array = P123_sparse_ptr_col_array[chn_3N];
        int     P123_dim       = P123_sparse_ptr_dim_array[chn_3N];

        /* CT and VC channel blocks */
        double** CT_RM_array = CT_RM_ptr_array[chn_3N];
        double** VC_CM_array = VC_CM_ptr_array[chn_3N];

        /* Convert P123_row_array from COO to CSR format */
        int* idx_row_array_csr = new int [block_dense_dim];
        coo_to_csr_format_converter(P123_row_array,
                                    idx_row_array_csr,
                                    P123_dim,
                                    block_dense_dim);
        P123_row_array = idx_row_array_csr;
        
        /* Create column of (C^T)(P)(VC) */
        double CPVC_col_array [block_dense_dim];


        /* Loop over columns of CPVC */
        printf("   - Calculating columns of CPVC ... \n");
        auto timestamp_calc_CPVC_start = std::chrono::system_clock::now();
        for (int idx_alpha_c=0; idx_alpha_c<Nalpha_block; idx_alpha_c++){
            for (int idx_q_c=0; idx_q_c<Nq_WP; idx_q_c++){
                for (int idx_p_c=0; idx_p_c<Np_WP; idx_p_c++){

                    /* Reset array */
                    for (int idx=0; idx<block_dense_dim; idx++){
                        CPVC_col_array[idx] = 0;
                    }

                    /* caluclate CPVC-column */
                    calculate_CPVC_col(CPVC_col_array,
                                       idx_alpha_c, idx_p_c, idx_q_c,
                                       Nalpha_block, Nq_WP, Np_WP,
                                       CT_RM_array,
                                       VC_CM_array,
                                       P123_val_array,
                                       P123_row_array,
                                       P123_col_array,
                                       P123_dim);
                }
            }
        }
        auto timestamp_calc_CPVC_end = std::chrono::system_clock::now();
		std::chrono::duration<double> time_calc_CPVC = timestamp_calc_CPVC_end - timestamp_calc_CPVC_start;
		printf("   - Done. Time used: %.6f\n", time_calc_CPVC.count());
    }
}

/* Create array of pointers to C^T matrices for product (C^T)PVC in row-major format */
void create_CT_row_maj_3N_pointer_array(double*** CT_RM_array,
                                        double*   C_WP_unco_array,
					                    double*   C_WP_coup_array,
                                        int       Np_WP,
                                        int       J_2N_max,
                                        int       N_chn_3N,
    	                                int*      chn_3N_idx_array,
					                    int       Nalpha,
                                        int*      L_2N_array,
                                        int*      S_2N_array,
                                        int*      J_2N_array,
                                        int*      T_2N_array){
    double* C_subarray  = NULL;
    double* CT_subarray = NULL;

    double* CT_unco_array = new double [Np_WP*Np_WP   * 2*(J_2N_max+1)];
    double* CT_coup_array = new double [Np_WP*Np_WP*4 *    J_2N_max   ];
    
    /* Copy and transpose all 2N-uncoupled C-arrays */
    for (int J_2N=0; J_2N<J_2N_max+1; J_2N++){
        for (int S_2N=0; S_2N<2; S_2N++){
            int idx_chn_unco       = 2*J_2N + S_2N;
			int idx_2N_mat_WP_unco = idx_chn_unco*Np_WP*Np_WP;
            
            C_subarray  = &C_WP_unco_array[idx_2N_mat_WP_unco];
            CT_subarray = &CT_unco_array  [idx_2N_mat_WP_unco];

            /* Copy content to avoid rewriting C-arrays */
            std::copy(C_subarray, C_subarray + Np_WP*Np_WP, CT_subarray);
            
            /* Transpose C to get C^T */
            simple_transpose_matrix_routine(CT_subarray, Np_WP);

        }
    }

    /* Copy and transpose all 2N-coupled C-arrays */
    for (int J_2N=1; J_2N<J_2N_max+1; J_2N++){
        int idx_chn_coup     = J_2N-1;
		int idx_2N_mat_WP_coup = idx_chn_coup*4*Np_WP*Np_WP;

        C_subarray  = &C_WP_coup_array[idx_2N_mat_WP_coup];
        CT_subarray = &CT_coup_array  [idx_2N_mat_WP_coup];
        
        /* Copy content to avoid rewriting C-arrays */
        std::copy(C_subarray, C_subarray + Np_WP*Np_WP, CT_subarray);
        
        /* Transpose C to get C^T */
        simple_transpose_matrix_routine(CT_subarray, 2*Np_WP);
        /* Restucture coupled matrix array into 4 seperate arrays */
        restructure_coupled_VC_product(CT_subarray, Np_WP);
    }
    
    /* 3N-channel loop */
    for (int chn_3N=0; chn_3N< N_chn_3N; chn_3N++){

        int idx_alpha_lower = chn_3N_idx_array[chn_3N];
        int idx_alpha_upper = chn_3N_idx_array[chn_3N+1];

        int Nalpha_block    = idx_alpha_upper - idx_alpha_lower;
        CT_RM_array[chn_3N] = new double* [Nalpha_block*Nalpha_block];

        /* Row state */
        for (int idx_alpha_r=idx_alpha_lower; idx_alpha_r<idx_alpha_upper; idx_alpha_r++){
            int L_r = L_2N_array[idx_alpha_r];
            int S_r = S_2N_array[idx_alpha_r];
            int J_r = J_2N_array[idx_alpha_r];
            int T_r = T_2N_array[idx_alpha_r];

            /* Column state */
            for (int idx_alpha_c=idx_alpha_lower; idx_alpha_c<idx_alpha_upper; idx_alpha_c++){
                int L_c = L_2N_array[idx_alpha_c];
                int S_c = S_2N_array[idx_alpha_c];
                int J_c = J_2N_array[idx_alpha_c];
                int T_c = T_2N_array[idx_alpha_c];

                /* Check if possible channel through interaction */
                if (T_r==T_c and J_r==J_c and S_r==S_c and abs(L_r-L_c)<=2){

                    /* Detemine if this is a coupled channel */
	                bool coupled_matrix = false;
	                if (L_r!=L_c or (L_r==L_c and L_r!=J_r and J_r!=0)){ // This counts 3P0 as uncoupled; used in matrix structure
	                	coupled_matrix  = true;
	                }

                    /* find which VC-product corresponds to the current coupling */
                    if (coupled_matrix){
                        int idx_chn_coup       = J_r-1;
	    	            int idx_2N_mat_WP_coup = idx_chn_coup*4*Np_WP*Np_WP;
                        if (L_r<L_c){       // L_r=J_r-1, L_c=J_r+1
                            CT_subarray = &CT_coup_array[idx_2N_mat_WP_coup + 1*Np_WP*Np_WP];
                        }
                        else if (L_r>L_c){  // L_r=J_r+1, L_c=J_r-1
                            CT_subarray = &CT_coup_array[idx_2N_mat_WP_coup + 2*Np_WP*Np_WP];
                        }
                        else if (L_r<J_c){  // L_r=J_r-1, L_c=J_r-1
                            CT_subarray = &CT_coup_array[idx_2N_mat_WP_coup + 0*Np_WP*Np_WP];
                        }
                        else{               // L_r=J_r+1, L_c=J_r+1
                            CT_subarray = &CT_coup_array[idx_2N_mat_WP_coup + 3*Np_WP*Np_WP];
                        }
                    }
                    else{
                        int idx_chn_unco       = 2*J_r + S_r;
	    		        int idx_2N_mat_WP_unco = idx_chn_unco*Np_WP*Np_WP;
                        CT_subarray = &CT_unco_array[idx_2N_mat_WP_unco];
                    }
                }
                else{
                    CT_subarray = NULL;
                }

                /* Unique index for two given states alpha */
                int idx_CT_RM = (idx_alpha_r-idx_alpha_lower)*Nalpha_block + (idx_alpha_c-idx_alpha_lower);
                (CT_RM_array[chn_3N])[idx_CT_RM] = CT_subarray;
            }
        }
    }
}

/* Restructures NN coupled matrix as 4 seperate matrices
 * !!! WARNING: COLUMN-MAJOR ALGORITHM !!! */
void restructure_coupled_VC_product(double* VC_product, int Np_WP){
    double* VC_product_temp = new double [4*Np_WP*Np_WP];

    for (int col_block=0; col_block<2; col_block++){
        for (int row_block=0; row_block<2; row_block++){

            int idx_block = col_block*2*Np_WP*Np_WP + row_block*Np_WP*Np_WP;

            for (int col=0; col<Np_WP; col++){
                for (int row=0; row<Np_WP; row++){
                    int prestructure_col_idx = (col + Np_WP*col_block)*2*Np_WP;
                    int prestructure_row_idx =  row + Np_WP*row_block;

                    VC_product_temp[idx_block + col*Np_WP+row] = VC_product[prestructure_col_idx + prestructure_row_idx];
                }
            }

        }
    }

    /* Copy from temp array */
    std::copy(VC_product_temp, VC_product_temp + 4*Np_WP*Np_WP, VC_product);

    /* Delete temporary array */
    delete [] VC_product_temp;
}

/* Create array of pointers to VC-product matrices for product (C^T)PVC in column-major format*/
void create_VC_col_maj_3N_pointer_array(double*** VC_CM_array,
                                        double*   C_WP_unco_array,
					                    double*   C_WP_coup_array,
					                    double*   V_WP_unco_array,
                                        double*   V_WP_coup_array,
					                    int       Np_WP,
                                        int       J_2N_max,
                                        int       N_chn_3N,
    	                                int*      chn_3N_idx_array,
					                    int       Nalpha,
                                        int*      L_2N_array,
                                        int*      S_2N_array,
                                        int*      J_2N_array,
                                        int*      T_2N_array){
    
    double* V_subarray = NULL;
    double* C_subarray = NULL;

    double* VC_unco_array = new double [Np_WP*Np_WP   * 2*(J_2N_max+1)];
    double* VC_coup_array = new double [Np_WP*Np_WP*4 *    J_2N_max   ];

    double* VC_product = NULL;

    /* Calculate all 2N-uncoupled VC-products and convert to column-major format */
    for (int J_2N=0; J_2N<J_2N_max+1; J_2N++){
        for (int S_2N=0; S_2N<2; S_2N++){
            int idx_chn_unco     = 2*J_2N + S_2N;
			int idx_2N_mat_WP_unco = idx_chn_unco*Np_WP*Np_WP;
            
            V_subarray = &V_WP_unco_array[idx_2N_mat_WP_unco];
            C_subarray = &C_WP_unco_array[idx_2N_mat_WP_unco];
            VC_product = &VC_unco_array  [idx_2N_mat_WP_unco];

            /* Multiply V and C using BLAS */
            dot_MM(V_subarray, C_subarray, VC_product, Np_WP, Np_WP, Np_WP);

            /* Transpose VC-product to get column-major format */
            simple_transpose_matrix_routine(VC_product, Np_WP);
        }
    }

    /* Calculate all 2N-coupled VC-products and convert to column-major format */
    for (int J_2N=1; J_2N<J_2N_max+1; J_2N++){
        int idx_chn_coup     = J_2N-1;
		int idx_2N_mat_WP_coup = idx_chn_coup*4*Np_WP*Np_WP;

        V_subarray = &V_WP_coup_array[idx_2N_mat_WP_coup];
        C_subarray = &C_WP_coup_array[idx_2N_mat_WP_coup];
        VC_product = &VC_coup_array  [idx_2N_mat_WP_coup];

        /* Multiply V and C using BLAS */
        dot_MM(V_subarray, C_subarray, VC_product, 2*Np_WP, 2*Np_WP, 2*Np_WP);
        /* Transpose VC-product to get column-major format */
        simple_transpose_matrix_routine(VC_product, 2*Np_WP);
        /* Restucture coupled matrix array into 4 seperate arrays */
        restructure_coupled_VC_product(VC_product, Np_WP);
    }

    /* 3N-channel loop */
    for (int chn_3N=0; chn_3N< N_chn_3N; chn_3N++){

        int idx_alpha_lower = chn_3N_idx_array[chn_3N];
        int idx_alpha_upper = chn_3N_idx_array[chn_3N+1];

        int Nalpha_block    = idx_alpha_upper - idx_alpha_lower;
        VC_CM_array[chn_3N] = new double* [Nalpha_block*Nalpha_block];

        /* Row state */
        for (int idx_alpha_r=idx_alpha_lower; idx_alpha_r<idx_alpha_upper; idx_alpha_r++){
            int L_r = L_2N_array[idx_alpha_r];
            int S_r = S_2N_array[idx_alpha_r];
            int J_r = J_2N_array[idx_alpha_r];
            int T_r = T_2N_array[idx_alpha_r];

            /* Column state */
            for (int idx_alpha_c=idx_alpha_lower; idx_alpha_c<idx_alpha_upper; idx_alpha_c++){
                int L_c = L_2N_array[idx_alpha_c];
                int S_c = S_2N_array[idx_alpha_c];
                int J_c = J_2N_array[idx_alpha_c];
                int T_c = T_2N_array[idx_alpha_c];

                /* Check if possible channel through interaction */
                if (T_r==T_c and J_r==J_c and S_r==S_c and abs(L_r-L_c)<=2){

                    /* Detemine if this is a coupled channel */
	                bool coupled_matrix = false;
	                if (L_r!=L_c or (L_r==L_c and L_r!=J_r and J_r!=0)){ // This counts 3P0 as uncoupled; used in matrix structure
	                	coupled_matrix  = true;
	                }

                    /* find which VC-product corresponds to the current coupling */
                    if (coupled_matrix){
                        int idx_chn_coup       = J_r-1;
	    	            int idx_2N_mat_WP_coup = idx_chn_coup*4*Np_WP*Np_WP;
                        if (L_r<L_c){       // L_r=J_r-1, L_c=J_r+1
                            VC_product = &VC_coup_array[idx_2N_mat_WP_coup + 1*Np_WP*Np_WP];
                        }
                        else if (L_r>L_c){  // L_r=J_r+1, L_c=J_r-1
                            VC_product = &VC_coup_array[idx_2N_mat_WP_coup + 2*Np_WP*Np_WP];
                        }
                        else if (L_r<J_c){  // L_r=J_r-1, L_c=J_r-1
                            VC_product = &VC_coup_array[idx_2N_mat_WP_coup + 0*Np_WP*Np_WP];
                        }
                        else{               // L_r=J_r+1, L_c=J_r+1
                            VC_product = &VC_coup_array[idx_2N_mat_WP_coup + 3*Np_WP*Np_WP];
                        }
                    }
                    else{
                        int idx_chn_unco       = 2*J_r + S_r;
	    		        int idx_2N_mat_WP_unco = idx_chn_unco*Np_WP*Np_WP;
                        VC_product = &VC_unco_array[idx_2N_mat_WP_unco];
                    }
                }
                else{
                    VC_product = NULL;
                }

                /* Unique index for two given states alpha */
                int idx_VC_CM = (idx_alpha_r-idx_alpha_lower)*Nalpha_block + (idx_alpha_c-idx_alpha_lower);
                (VC_CM_array[chn_3N])[idx_VC_CM] = VC_product;
            }
        }
    }
}
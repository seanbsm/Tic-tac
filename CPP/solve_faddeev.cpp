
#include "solve_faddeev.h"

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


void iterate_faddeev_equations(cdouble* U_array,
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
    
}

/* Solves the Faddeev equations
 * U = P*V + P*V*G*U
 * on the form L*U = R, where L and R are the left-
 * and right-handed sides of the equations, given by 
 * L = 1 - P*V*G
 * R = P*V
 * Since G is expressed in an SWP basis, we also must include the basis-transormation matrices C */
void direct_solve_faddeev_equations(cdouble* U_array,
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
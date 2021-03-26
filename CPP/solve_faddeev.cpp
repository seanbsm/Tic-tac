
#include "solve_faddeev.h"

void calculate_CPVC_product_row(){
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
						double VC_element = VC_subarray_col[idx_p_j];
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
							double PVC_element = PVC_col[idx_PVC];

							/* I'm not sure if this is the fastest ordering of the loops */
							double CT_element  = CT_subarray_row[idx_p_i];

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

void CPVC_col_brute_force(double*  col_array,
						  int      idx_alpha_c, int idx_p_c, int idx_q_c,
						  int      Nalpha,      int Nq_WP,   int Np_WP,
						  double** CT_RM_array,
						  double** VC_CM_array,
						  double*  P123_val_array,
						  int*     P123_row_array,
						  int*     P123_col_array,
						  int      P123_dim){
	double* CT_ptr = NULL;
	double* VC_ptr = NULL;

	bool print_content = false;

	/* Index: r */
	for (int idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
		for (int idx_q_r=0; idx_q_r<Nq_WP; idx_q_r++){
			for (int idx_p_r=0; idx_p_r<Np_WP; idx_p_r++){
				if (print_content){std::cout << "Index i " << idx_alpha_r << " " << idx_q_r << " " << idx_p_r << std::endl;}
				double sum = 0;
				/* Index: i */
				for (int idx_alpha_i=0; idx_alpha_i<Nalpha; idx_alpha_i++){
					CT_ptr = CT_RM_array[idx_alpha_r*Nalpha + idx_alpha_i];
					if (CT_ptr!=NULL){
						for (int idx_q_i=0; idx_q_i<Nq_WP; idx_q_i++){
							if (idx_q_r == idx_q_i){
								for (int idx_p_i=0; idx_p_i<Np_WP; idx_p_i++){
									if (print_content){std::cout << " - Index i " << idx_alpha_i << " " << idx_q_i << " " << idx_p_i << std::endl;}
									/* Index: j */
									for (int idx_alpha_j=0; idx_alpha_j<Nalpha; idx_alpha_j++){
										VC_ptr = VC_CM_array[idx_alpha_c*Nalpha + idx_alpha_j];
										if (VC_ptr!=NULL){
											for (int idx_q_j=0; idx_q_j<Nq_WP; idx_q_j++){
												if (idx_q_j == idx_q_c){
													for (int idx_p_j=0; idx_p_j<Np_WP; idx_p_j++){
														if (print_content){std::cout << "  - Index j " << idx_alpha_j << " " << idx_q_j << " " << idx_p_j << std::endl;}

														int idx_P123_row = idx_alpha_i*Nq_WP*Np_WP + idx_q_i*Np_WP + idx_p_i;
														int idx_P123_col = idx_alpha_j*Nq_WP*Np_WP + idx_q_j*Np_WP + idx_p_j;

														int idx_P123_col_lower = P123_row_array[idx_P123_row];
														int idx_P123_col_upper = P123_row_array[idx_P123_row+1];
														bool idx_found = false;
														int idx_P123_val = 0;
														for (int nnz_idx=idx_P123_col_lower; nnz_idx<idx_P123_col_upper; nnz_idx++){
															if (P123_col_array[nnz_idx] == idx_P123_col){
																idx_found = true;
																idx_P123_val = nnz_idx;
															}
														}

														if (idx_found){
															double CT_element = CT_ptr[idx_p_r*Np_WP + idx_p_i];
															double P_element  = P123_val_array[idx_P123_val];
															double VC_element = VC_ptr[idx_p_c*Np_WP + idx_p_j];

															sum += CT_element * P_element * VC_element;
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}

				int idx_CPVC = idx_alpha_r*Nq_WP*Np_WP + idx_q_r*Np_WP + idx_p_r;
				col_array[idx_CPVC] = sum;
			}
		}
	}
}

void CPVC_col_calc_test(int      Nalpha,
						int 	 Nq_WP,
						int 	 Np_WP,
						double** CT_RM_array,
						double** VC_CM_array,
						double*  P123_sparse_val_array,
						int*     P123_sparse_row_array,
						int*     P123_sparse_col_array,
						int      P123_sparse_dim){

	bool print_content = false;

	int dense_dim = Nalpha * Nq_WP * Np_WP;

	/* Create column of (C^T)(P)(VC) */
	double CPVC_col_array    [dense_dim];
	double CPVC_col_array_BF [dense_dim];

	for (int idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
		for (int idx_q_c=0; idx_q_c<Nq_WP; idx_q_c++){
			for (int idx_p_c=0; idx_p_c<Np_WP; idx_p_c++){
	
				/* Reset array */
				for (int idx=0; idx<dense_dim; idx++){
					CPVC_col_array[idx]    = 0;
					CPVC_col_array_BF[idx] = 0;
				}
	
				/* Calculate CPVC-column */
				calculate_CPVC_col(CPVC_col_array,
								   idx_alpha_c, idx_p_c, idx_q_c,
								   Nalpha, Nq_WP, Np_WP,
								   CT_RM_array,
								   VC_CM_array,
								   P123_sparse_val_array,
								   P123_sparse_row_array,
								   P123_sparse_col_array,
								   P123_sparse_dim);
				
				CPVC_col_brute_force(CPVC_col_array_BF,
									 idx_alpha_c, idx_p_c, idx_q_c,
									 Nalpha, Nq_WP, Np_WP,
									 CT_RM_array,
									 VC_CM_array,
									 P123_sparse_val_array,
									 P123_sparse_row_array,
									 P123_sparse_col_array,
									 P123_sparse_dim);
				
				for (int idx=0; idx<dense_dim; idx++){
					double diff = abs(CPVC_col_array[idx] - CPVC_col_array_BF[idx]);
					if ( diff > 1e-14 ){
						printf("Element %d had a discrepency: %.16f vs. %.16f \n", idx, CPVC_col_array[idx], CPVC_col_array_BF[idx]);
						raise_error("CPVC benchmarking failed");
					}
				}
			}
		}
	}
}

cdouble pade_approximant(cdouble* a_coeff_array, int N, int M, cdouble z){

	/* a_coeff_array must have length N+M+1 */
	cdouble P_array [(M+1)*(M+1)];
	cdouble Q_array [(M+1)*(M+1)];

	for (int row_idx=0; row_idx<M; row_idx++){
		for (int col_idx=0; col_idx<M+1; col_idx++){
			P_array[row_idx*(M+1) + col_idx] = a_coeff_array[M-N+1 + row_idx + col_idx];
			Q_array[row_idx*(M+1) + col_idx] = a_coeff_array[M-N+1 + row_idx + col_idx];
		}
	}

	for (int col_idx=0; col_idx<M+1; col_idx++){
		Q_array[M*(M+1) + col_idx] = std::pow(z, M-col_idx);
		P_array[M*(M+1) + col_idx] = 0;
		for (int j=M-col_idx; j<N+1; j++){
			P_array[M*(M+1) + col_idx] += a_coeff_array[col_idx - M] * std::pow(z, j);
		}
	}

	cdouble P_det = determinant(P_array, M+1);
	cdouble Q_det = determinant(Q_array, M+1);

	return P_det/Q_det;
}

void pade_method_solve(cdouble*  U_array,
					   cdouble*  G_array,
					   int		 num_on_shell_indices,
					   int*      on_shell_idx_array,
					   int       Nalpha,
					   int 	     Nq_WP,
					   int 	     Np_WP,
					   double**  CT_RM_array,
					   double**  VC_CM_array,
					   double*   P123_sparse_val_array,
					   int*      P123_sparse_row_array,
					   int*      P123_sparse_col_array,
					   int       P123_sparse_dim){
	
	/* Upper limit on polynomial approximation of Faddeev eq. */
	int N_pade = 10;
	int M_pade = 10;

	/* Coefficients for calculating Pade approximant */
	cdouble* a_coeff_array = new cdouble [N_pade + M_pade + 1];

	/* Dense dimension of 3N-channel */
	int dense_dim = Nalpha * Nq_WP * Np_WP;

	/* Allocate column-array for (C^T)(P)(VC) */
	double CPVC_col_array [dense_dim];

	/* Loop over number of Pade-terms we use */
	for (int n=0; n<N_pade+M_pade; n++){
		/* Iterate through columns of A */
		for (int idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
			for (int idx_q_c=0; idx_q_c<Nq_WP; idx_q_c++){
				for (int idx_p_c=0; idx_p_c<Np_WP; idx_p_c++){
				
					/* Reset CPVC-column array */
					for (int row_idx=0; row_idx<dense_dim; row_idx++){
						CPVC_col_array[row_idx] = 0;
					}
	
					/* Calculate CPVC-column */
					calculate_CPVC_col(CPVC_col_array,
									   idx_alpha_c, idx_p_c, idx_q_c,
									   Nalpha, Nq_WP, Np_WP,
									   CT_RM_array,
									   VC_CM_array,
									   P123_sparse_val_array,
									   P123_sparse_row_array,
									   P123_sparse_col_array,
									   P123_sparse_dim);
					
					/* Multiply CPVC with resolvent, raised to the power n given by outer for-loop */
					cdouble G_pow_n = std::pow(G_array[idx_alpha_c*Nq_WP*Np_WP + idx_q_c*Np_WP + idx_p_c], n);

					/* Calculate all a-coefficients for calculated CPVC-column */
					//for (int idx=0; idx<num_on_shell_indices; idx++){
					//}
				}
			}
		}
	}
}

/* Solves Faddeev on the form
 * (1-AG)U = A
 * where A = C^T PVC */
void direct_sparse_solve(cdouble*  U_array,
						 cdouble*  G_array,
						 int       idx_on_shell,
						 int       Nalpha,
						 int 	   Nq_WP,
						 int 	   Np_WP,
						 double**  CT_RM_array,
						 double**  VC_CM_array,
						 double*   P123_sparse_val_array,
						 int*      P123_sparse_row_array,
						 int*      P123_sparse_col_array,
						 int       P123_sparse_dim){
	
	/* Dense dimension of 3N-channel */
	int dense_dim = Nalpha * Nq_WP * Np_WP;

	/* Sparse dimension and step-length for incrementing sparse array size */
	int		A_sparse_dim	   = 0;
	int 	sparse_step_length = 0;
	if (dense_dim>1000){
		sparse_step_length = dense_dim/1000;
	}
	else{
		sparse_step_length = dense_dim;
	}
	int current_array_dim  = sparse_step_length;

	/* Dynamically sized sparse-storage COO (CM) array format for A-matrix */
	double* A_sparse_val_array = new double [sparse_step_length];
	int* 	A_sparse_row_array = new int    [sparse_step_length];
	int* 	A_sparse_col_array = new int    [sparse_step_length];

	/* Allocate column-array for (C^T)(P)(VC) */
	double 				 CPVC_col_array 	  [dense_dim];

	std::complex<double> A_on_shell_col_array [dense_dim];
	std::complex<double>* A_dense_array = new std::complex<double> [dense_dim*dense_dim];
	//double A_on_shell_col_array [dense_dim];
	//double* A_dense_array	   = new double [dense_dim*dense_dim];
	
	for (int idx=0; idx<dense_dim*dense_dim; idx++){
		A_dense_array[idx] = 0;
	}

	/* Calculation of columns of A */
	for (int idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
		for (int idx_q_c=0; idx_q_c<Nq_WP; idx_q_c++){
			for (int idx_p_c=0; idx_p_c<Np_WP; idx_p_c++){

				/* Reset CPVC-column array */
				for (int row_idx=0; row_idx<dense_dim; row_idx++){
					CPVC_col_array[row_idx]    = 0;
				}

				/* Calculate CPVC-column */
				calculate_CPVC_col(CPVC_col_array,
								   idx_alpha_c, idx_p_c, idx_q_c,
								   Nalpha, Nq_WP, Np_WP,
								   CT_RM_array,
								   VC_CM_array,
								   P123_sparse_val_array,
								   P123_sparse_row_array,
								   P123_sparse_col_array,
								   P123_sparse_dim);
				
				/* Append on-shell column to right-hand side vector of Faddeev eq. (and convert type to complex) */
				int col_idx = idx_alpha_c*Np_WP*Nq_WP + idx_q_c*Np_WP + idx_p_c;
				if (col_idx == idx_on_shell){
					for (int row_idx=0; row_idx<dense_dim; row_idx++){
						A_on_shell_col_array[row_idx] = CPVC_col_array[row_idx];
					}
				}
				
				/* Loop through rows of CPVC-column and append to A if non-zero */
				for (int row_idx=0; row_idx<dense_dim; row_idx++){

					double element = CPVC_col_array[row_idx];
					/* Add element if non-zero or if this a diagonal element (needed for "1-AG" in Faddeev eq.) */
					if (element!=0 || col_idx==row_idx){

						/* Append to sparse value and index arrays */
						A_sparse_val_array[A_sparse_dim] = element;
						A_sparse_row_array[A_sparse_dim] = row_idx;
						A_sparse_col_array[A_sparse_dim] = col_idx;

						A_dense_array[row_idx*dense_dim + col_idx] = element;

						/* Increment sparse dimension (num of non-zero elements) */
						A_sparse_dim += 1;

						/* If the dimension goes over the array dimension we increase the array size
						 * via a copy-paste-type routine, and increment the current array dimension */
						if ( A_sparse_dim>=current_array_dim ){  // This should occur a small amount of the time
							increase_sparse_array_size(&A_sparse_val_array, current_array_dim, sparse_step_length);
							increase_sparse_array_size(&A_sparse_row_array, current_array_dim, sparse_step_length);
							increase_sparse_array_size(&A_sparse_col_array, current_array_dim, sparse_step_length);

							/* Increment sparse-array dimension */
							current_array_dim += sparse_step_length;
						}
					}
				}
			}
		}
	}
	/* Contract arrays to minimal size (number of non-zero elements) */
	reduce_sparse_array_size(&A_sparse_val_array, current_array_dim, A_sparse_dim);
	reduce_sparse_array_size(&A_sparse_row_array, current_array_dim, A_sparse_dim);
	reduce_sparse_array_size(&A_sparse_col_array, current_array_dim, A_sparse_dim);

	/* Conversion from column-major COO to row-major COO array format for A-matrix 
	 * !!! WARNING: THIS ROUTINE REWRITES INPUT ARRAYS TO NEW FORMAT, OLD FORMAT IS DELETED !!! */
	coo_col_major_to_coo_row_major_converter(&A_sparse_val_array,
											 &A_sparse_row_array,
											 &A_sparse_col_array,
											 A_sparse_dim,
											 dense_dim);
	/* Conversion from COO to CSR array format for A-matrix */
	int* A_idx_row_array_csr = new int [dense_dim+1];
	coo_to_csr_format_converter(A_sparse_row_array,
                                A_idx_row_array_csr,
                                A_sparse_dim,
                                dense_dim);

	/* Test sparse format */
	if (false){
		for (int row_idx=0; row_idx<dense_dim; row_idx++){
			int nnz_idx_lower = A_idx_row_array_csr[row_idx];
			int nnz_idx_upper = A_idx_row_array_csr[row_idx+1];
			for (int nnz_idx=nnz_idx_lower; nnz_idx<nnz_idx_upper; nnz_idx++){
				int col_idx = A_sparse_col_array[nnz_idx];

				//double val_sparse = A_sparse_val_array[nnz_idx];
				//double val_dense  = A_dense_array[row_idx*dense_dim + col_idx];
				std::complex<double> val_sparse = A_sparse_val_array[nnz_idx];
				std::complex<double> val_dense  = A_dense_array[row_idx*dense_dim + col_idx];
				if (val_sparse!=val_dense){
					std::cout << "Mismatch. Sparse val: " << val_sparse.real() << " " << val_sparse.imag() << std::endl;
					std::cout << "Mismatch. Dense val:  " << val_dense.real() << " " << val_dense.imag() << std::endl;
					raise_error("Sparse format conversion failed");
				}
			}
		}
	}

	/* Identity "I" minus "AG"-product and conversion to complex type */
	std::complex<double>* IAG_sparse_val_array_cmplx = new std::complex<double> [A_sparse_dim];
	std::complex<double>* IAG_dense_array_cmplx = new std::complex<double> [dense_dim*dense_dim];
	//double* IAG_sparse_val_array_cmplx = new double [A_sparse_dim];
	//double* IAG_dense_array_cmplx	   = new double [dense_dim*dense_dim];
	for (int idx=0; idx<dense_dim*dense_dim; idx++){
		IAG_dense_array_cmplx[idx] = 0;
	}
	for (int row_idx=0; row_idx<dense_dim; row_idx++){

		int nnz_idx_lower = A_idx_row_array_csr[row_idx];
		int nnz_idx_upper = A_idx_row_array_csr[row_idx+1];

		for (int nnz_idx=nnz_idx_lower; nnz_idx<nnz_idx_upper; nnz_idx++){

			int col_idx = A_sparse_col_array[nnz_idx];

			std::complex<double> G_val = G_array[col_idx];

			/* AG-multiplication and minus-sign */
			IAG_sparse_val_array_cmplx[nnz_idx] = -A_sparse_val_array[nnz_idx] * G_val;
			//IAG_sparse_val_array_cmplx[nnz_idx] = -A_sparse_val_array[nnz_idx] * G_val.real();

			/* Add identity matrix */
			if (row_idx==col_idx){
				IAG_sparse_val_array_cmplx[nnz_idx] += 1;
			}

			/* Dense test case */
			IAG_dense_array_cmplx[row_idx*dense_dim + col_idx] = IAG_sparse_val_array_cmplx[nnz_idx];
		}
	}

	MKL_INT rows [dense_dim+1];
	for (int row_idx=0; row_idx<dense_dim+1; row_idx++){
		rows[row_idx] = A_idx_row_array_csr[row_idx];
	}
	MKL_INT cols [A_sparse_dim];
	for (int col_idx=0; col_idx<A_sparse_dim; col_idx++){
		cols[col_idx] = A_sparse_col_array[col_idx];
	}

	MKL_INT sparse_dim_MKL = A_sparse_dim;
	MKL_INT dense_dim_MKL = dense_dim;

	//double sols_U_col [dense_dim];
	std::complex<double> sols_U_col [dense_dim];

	/* Solve Faddeev using MKL DSS-routines */
	auto timestamp_1 = std::chrono::system_clock::now();
	solve_MM_sparse(IAG_sparse_val_array_cmplx,
                 	rows,//(MKL_INT*) A_idx_row_array_csr,
                 	cols,//(MKL_INT*) A_sparse_col_array,
                 	sparse_dim_MKL,//(MKL_INT) A_sparse_dim,
                 	A_on_shell_col_array,
                 	dense_dim_MKL,//(MKL_INT) dense_dim,
					sols_U_col);
	auto timestamp_2 = std::chrono::system_clock::now();
	solve_MM(IAG_dense_array_cmplx, A_dense_array, dense_dim);
	auto timestamp_3 = std::chrono::system_clock::now();
	std::chrono::duration<double> time_sparse = timestamp_2 - timestamp_1;
	std::chrono::duration<double> time_dense  = timestamp_3 - timestamp_2;
	printf("   - Time used sparse: %.6f\n", time_sparse.count());
	printf("   - Time used dense:  %.6f\n", time_dense.count());

	for (int row=0; row<dense_dim; row++){
		//double val_sparse = sols_U_col[row];
		//double val_dense  = IAG_dense_array_cmplx[row*dense_dim + idx_on_shell];
		std::complex<double> val_sparse = sols_U_col[row];
		std::complex<double> val_dense  = IAG_dense_array_cmplx[row*dense_dim + idx_on_shell];
		//if (val_sparse!=val_dense and std::isnan(val_dense.real())!=true and std::isnan(val_sparse.real())!=true
		//and std::isnan(val_dense.imag())!=true and std::isnan(val_sparse.imag())!=true){
			std::cout << "Mismatch. Sparse val: " << val_sparse.real() << " " << val_sparse.imag() << std::endl;
			std::cout << "Mismatch. Dense val:  " << val_dense.real() << " " << val_dense.imag() << std::endl;
		//	raise_error("Exit");
		//}
	}

	/* Delete temporary arrays */
	delete [] A_sparse_val_array;
	delete [] A_sparse_row_array;
	delete [] A_sparse_col_array;
	delete [] IAG_sparse_val_array_cmplx;
}

void solve_faddeev_equations(cdouble*  U_array,
							 cdouble*  G_array,
							 double*   P123_sparse_val_array,
							 int*      P123_sparse_row_array,
							 int*      P123_sparse_col_array,
							 int       P123_sparse_dim,
							 double*   C_WP_unco_array,
							 double*   C_WP_coup_array,
							 double*   V_WP_unco_array,
							 double*   V_WP_coup_array,
							 int       idx_on_shell,
							 int       J_2N_max,
							 int       Nq_WP,
							 int       Np_WP,
							 int       Nalpha,
							 int*      L_2N_array,
							 int*      S_2N_array,
							 int*      J_2N_array,
							 int*      T_2N_array,
							 int*      L_1N_array, 
							 int*      two_J_1N_array){
	
	bool test_CPVC_col_routine  = false;
	bool use_DSS_solver_routine = false;

	/* Create C^T-product pointer-arrays in row-major format */
	double** CT_RM_array = new double* [Nalpha*Nalpha];
	create_CT_row_maj_3N_pointer_array(CT_RM_array,
									   C_WP_unco_array,
									   C_WP_coup_array,
									   Np_WP,
									   J_2N_max,
									   Nalpha,
									   L_2N_array,
									   S_2N_array,
									   J_2N_array,
									   T_2N_array);

	/* Create VC-product pointer-arrays in column-major format */
	double** VC_CM_array  = new double* [Nalpha*Nalpha];
	create_VC_col_maj_3N_pointer_array(VC_CM_array,
									   C_WP_unco_array,
									   C_WP_coup_array,
									   V_WP_unco_array,
									   V_WP_coup_array,
									   Np_WP,
									   J_2N_max,
									   Nalpha,
									   L_2N_array,
									   S_2N_array,
									   J_2N_array,
									   T_2N_array);
	
	/* Convert P123_row_array from COO to CSR format */
	int dense_dim = Nalpha * Nq_WP * Np_WP;
	int* idx_row_array_csr = new int [dense_dim];
	coo_to_csr_format_converter(P123_sparse_row_array,
								idx_row_array_csr,
								P123_sparse_dim,
								dense_dim);
	P123_sparse_row_array = idx_row_array_csr;
		
	/* Test optimized routine for CPVC columns */
	if (test_CPVC_col_routine){
		CPVC_col_calc_test(Nalpha,
						   Nq_WP,
						   Np_WP,
						   CT_RM_array,
						   VC_CM_array,
						   P123_sparse_val_array,
						   P123_sparse_row_array,
						   P123_sparse_col_array,
						   P123_sparse_dim);
	}
	
	printf("   - Solving Faddeev equation ... \n");
	auto timestamp_calc_CPVC_start = std::chrono::system_clock::now();
	
	direct_sparse_solve(U_array,
						G_array,
						idx_on_shell,
						Nalpha,
						Nq_WP,
						Np_WP,
						CT_RM_array,
						VC_CM_array,
						P123_sparse_val_array,
						P123_sparse_row_array,
						P123_sparse_col_array,
						P123_sparse_dim);

	auto timestamp_calc_CPVC_end = std::chrono::system_clock::now();
	std::chrono::duration<double> time_calc_CPVC = timestamp_calc_CPVC_end - timestamp_calc_CPVC_start;
	printf("   - Done. Time used: %.6f\n", time_calc_CPVC.count());
}

/* Create array of pointers to C^T matrices for product (C^T)PVC in row-major format */
void create_CT_row_maj_3N_pointer_array(double** CT_RM_array,
										double*  C_WP_unco_array,
										double*  C_WP_coup_array,
										int      Np_WP,
										int      J_2N_max,
										int      Nalpha,
										int*     L_2N_array,
										int*     S_2N_array,
										int*     J_2N_array,
										int*     T_2N_array){

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
			int idx_CT_RM = idx_alpha_r*Nalpha + idx_alpha_c;
			CT_RM_array[idx_CT_RM] = CT_subarray;
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
void create_VC_col_maj_3N_pointer_array(double** VC_CM_array,
										double*  C_WP_unco_array,
										double*  C_WP_coup_array,
										double*  V_WP_unco_array,
										double*  V_WP_coup_array,
										int      Np_WP,
										int      J_2N_max,
										int      Nalpha,
										int*     L_2N_array,
										int*     S_2N_array,
										int*     J_2N_array,
										int*     T_2N_array){

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
			int idx_VC_CM = idx_alpha_r*Nalpha + idx_alpha_c;
			VC_CM_array[idx_VC_CM] = VC_product;
		}
	}
}
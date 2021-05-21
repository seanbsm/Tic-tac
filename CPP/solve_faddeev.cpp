
#include "solve_faddeev.h"

void calculate_PVC_col(double*  col_array,
					   size_t   idx_alpha_c, size_t idx_p_c, size_t idx_q_c,
					   size_t   Nalpha,      size_t Nq_WP,   size_t Np_WP,
					   double** VC_CM_array,
					   double*  P123_val_array,
					   int*  	P123_row_array,
					   size_t*  P123_col_array,
					   size_t   P123_dim){

	double* VC_subarray     = NULL;

	size_t dense_dim = Nalpha*Np_WP*Nq_WP;

	for (size_t idx_alpha_j=0; idx_alpha_j<Nalpha; idx_alpha_j++){
		VC_subarray = VC_CM_array[idx_alpha_c*Nalpha + idx_alpha_j];

		/* Only do inner-product if VC is not zero due to conservation laws */
		if (VC_subarray!=NULL){
			for (size_t idx_p_j=0; idx_p_j<Np_WP; idx_p_j++){

				/* Access VC element */
				double VC_element = VC_subarray[idx_p_c*Np_WP + idx_p_j];

				/* Inner-product index */
				size_t idx_j =  idx_alpha_j*Np_WP*Nq_WP + idx_q_c*Np_WP + idx_p_j;

				/* CSC-format indexing */
				size_t idx_i_lower = P123_col_array[idx_j    ];
				size_t idx_i_upper = P123_col_array[idx_j + 1];

				/* Loop through rows of column we're calculating, and append */
				for (size_t idx_i=idx_i_lower; idx_i<idx_i_upper; idx_i++){

					/* Access arrays, this whole function is written to minimize these two calls */
					/* NOTE THAT P = P123 + P132 = 2*P123 */
					double P_element  = 2*P123_val_array[idx_i];

					size_t idx_row = P123_row_array[idx_i];
					/* Write inner product to col_array of PVC-product */
					col_array[idx_row] += P_element * VC_element;
				}
			}
		}
	}
}

void calculate_CPVC_col(double*  col_array,
						int* 	 row_to_nnz_array, 
						int* 	 nnz_to_row_array,
						size_t&  num_nnz,
						size_t   idx_alpha_c, size_t idx_p_c, size_t idx_q_c,
						size_t   Nalpha,      size_t Nq_WP,   size_t Np_WP,
						double** CT_RM_array,
						double** VC_CM_array,
						double*  P123_val_array,
						int*     P123_row_array,
						size_t*  P123_col_array,
						size_t   P123_dim){
	
	/* Generate PVC-column */
	double* PVC_col = new double [Nalpha*Nq_WP*Np_WP];
	/* Ensure PVC_col contains only zeroes */
	for (size_t idx=0; idx<Nalpha*Nq_WP*Np_WP; idx++){
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

	/* THOUGHT:
	 * MOVE ALPHA_I OUTWARDS AND GO BACK TO DIRECT APPEND TO COL_ARRAY.
	 * BUT USE NON-ZERO IF-TESTING. */

	/* TOUGHT (CONTRADICTORY TO DIRECT APPEND):
	 * Use dense mat-vec multiplication for sub-blocks  (Can let q be columns of right-vectors? Appealing use of MM-multiplication)
	 * Use dense vec-vec (or vec-mat-vec?) multiplication for A_An */

	//#pragma omp parallel
	//{
	/* Generate (C^T x PVC)-column */
	double* CT_subarray     = NULL;
	double* CT_subarray_row = NULL;
	double* PVC_subcol 		= NULL;
	//#pragma omp for
	/* Loop over rows of col_array */
	for (size_t idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
		/* Beginning of inner-product loops (index "i") */
		for (size_t idx_alpha_i=0; idx_alpha_i<Nalpha; idx_alpha_i++){
			size_t idx_CT_2N_block = idx_alpha_r*Nalpha + idx_alpha_i;
			CT_subarray = CT_RM_array[idx_CT_2N_block];

			/* Only do inner-product if CT is not zero due to conservation laws */
			if (CT_subarray!=NULL){
				PVC_subcol = &PVC_col[idx_alpha_i*Nq_WP*Np_WP];
				for (size_t idx_q_r=0; idx_q_r<Nq_WP; idx_q_r++){
					for (size_t idx_p_i=0; idx_p_i<Np_WP; idx_p_i++){
						//size_t idx_PVC     = idx_alpha_i*Nq_WP*Np_WP + idx_q_r*Np_WP + idx_p_i;
						double PVC_element = PVC_subcol[idx_q_r*Np_WP + idx_p_i];
						if (PVC_element!=0){
							for (size_t idx_p_r=0; idx_p_r<Np_WP; idx_p_r++){

								/* I'm not sure if this is the fastest ordering of the loops */
								double CT_element  = CT_subarray[idx_p_r*Np_WP + idx_p_i];

								double prod = CT_element * PVC_element;
								if (prod!=0){
									size_t idx_CPVC = idx_alpha_r*Nq_WP*Np_WP + idx_q_r*Np_WP + idx_p_r;
									//col_array[idx_CPVC] += prod;

									int nnz_idx = row_to_nnz_array[idx_CPVC];
									if (nnz_idx==-1){
										row_to_nnz_array[idx_CPVC] = num_nnz;
										nnz_to_row_array[num_nnz]  = idx_CPVC;
										nnz_idx = num_nnz;
										num_nnz += 1;
									}
									col_array[nnz_idx] += prod;
								}
							}
						}
					}
				}
			}
		}
	}
	//}
	delete [] PVC_col;
}

void calculate_all_CPVC_rows(cdouble* row_arrays,
							 int*	  q_com_idx_array,	  size_t num_q_com,
					   		 int*     deuteron_idx_array, size_t num_deuteron_states,
							 size_t   Nalpha, 
							 size_t   Nq_WP,
							 size_t   Np_WP,
							 double** CT_RM_array,
							 double** VC_CM_array,
							 double*  P123_val_array,
							 int*     P123_row_array,
							 size_t*  P123_col_array,
							 size_t   P123_dim){
	
	size_t dense_dim = Nalpha*Nq_WP*Np_WP;

	/* Loop over cols of row */
	#pragma omp parallel
	{
	/* Generate PVC-column */
	double* PVC_col = new double [dense_dim];
	/* Generate (C^T x PVC)-column */
	double* CT_subarray     = NULL;
	double* CT_subarray_row = NULL;
	#pragma omp for
	for (size_t idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
		for (size_t idx_q_c=0; idx_q_c<Nq_WP; idx_q_c++){
			for (size_t idx_p_c=0; idx_p_c<Np_WP; idx_p_c++){

				/* Ensure PVC_col contains only zeroes */
				for (size_t idx=0; idx<dense_dim; idx++){
					PVC_col[idx] = 0;
				}
				
				/* Calculate PVC-column for alpha_i, p_i, q_r */
				calculate_PVC_col(PVC_col,
								  idx_alpha_c, idx_p_c, idx_q_c,
								  Nalpha,      Nq_WP,   Np_WP,
								  VC_CM_array,
								  P123_val_array,
								  P123_row_array,
								  P123_col_array,
								  P123_dim);

				/* Re-use PVC-column in all relevant calculations */
				for (size_t i=0; i<num_deuteron_states; i++){
					for (size_t j=0; j<num_q_com; j++){
						/* Nucleon-deuteron on-shell (NDOS) indices
						 * (deuteron bound-state p-index is alwasy 0 due to eigenvalue ordering in SWP construction) */
						size_t idx_alpha_r = deuteron_idx_array[i];
						size_t idx_p_r     = 0;
						size_t idx_q_r     = q_com_idx_array[j];

						size_t idx_NDOS = i*num_q_com + j;;

						double inner_product_CPVC = 0;
						/* Beginning of inner-product loops (index "i") */
						for (size_t idx_alpha_i=0; idx_alpha_i<Nalpha; idx_alpha_i++){
							
							CT_subarray = CT_RM_array[idx_alpha_r*Nalpha + idx_alpha_i];
		
							/* Only do inner-product if CT is not zero due to conservation laws */
							if (CT_subarray!=NULL){
								CT_subarray_row = &CT_subarray[idx_p_r*Np_WP];
		
								for (size_t idx_p_i=0; idx_p_i<Np_WP; idx_p_i++){
								
									size_t idx_PVC     = idx_alpha_i*Nq_WP*Np_WP + idx_q_r*Np_WP + idx_p_i;
									double PVC_element = PVC_col[idx_PVC];
		
									/* I'm not sure if this is the fastest ordering of the loops */
									double CT_element  = CT_subarray_row[idx_p_i];
		
									inner_product_CPVC += CT_element * PVC_element;
								}
							}
						}
				
						size_t idx_CPVC = idx_alpha_c*Nq_WP*Np_WP + idx_q_c*Np_WP + idx_p_c;
						row_arrays[idx_NDOS*dense_dim + idx_CPVC] = inner_product_CPVC;
					}
				}
			}
		}
	}
	delete [] PVC_col; 
	}
}

cdouble pade_approximant(cdouble* a_coeff_array, size_t N, size_t M, cdouble z){

	/* a_coeff_array must have length N+M+1 */
	cdouble P_array [(M+1)*(M+1)];
	cdouble Q_array [(M+1)*(M+1)];

	for (size_t row_idx=0; row_idx<M; row_idx++){
		for (size_t col_idx=0; col_idx<M+1; col_idx++){
			P_array[row_idx*(M+1) + col_idx] = a_coeff_array[N-M+1 + row_idx + col_idx];
			Q_array[row_idx*(M+1) + col_idx] = a_coeff_array[N-M+1 + row_idx + col_idx];
		}
	}

	for (size_t col_idx=0; col_idx<M+1; col_idx++){
		Q_array[M*(M+1) + col_idx] = std::pow(z, M-col_idx);
		P_array[M*(M+1) + col_idx] = 0;
		for (size_t j=M-col_idx; j<N+1; j++){
			P_array[M*(M+1) + col_idx] += a_coeff_array[j - (M-col_idx)] * std::pow(z, j);
		}
	}

	cdouble P_det = determinant(P_array, M+1);
	cdouble Q_det = determinant(Q_array, M+1);

	return P_det/Q_det;
}

#include <fstream>
#include <iomanip>
void store_array(cdouble* array, size_t array_length, std::string filename){
	/* Open file*/
	std::ofstream result_file;
	result_file.open(filename);
	
	/* Fixes formatting of stored numbers */
	result_file << std::fixed
				<< std::showpos
				//~ << std::right
                //~ << std::setw(14)
				<< std::setprecision(8);
	
	/* Append array-values */
	for (size_t i=0; i<array_length; i++){
		/* Append vector element */
		//result_file << array[i] << "\n";
        result_file << "(" << array[i].real() << array[i].imag() << "j)\n";
	}
	
	/* Close writing session */
	result_file << std::endl;
	
	/* Close files */
	result_file.close();
}

/* Solves the Faddeev equations
 * U = P*V + P*V*G*U
 * on the form L*U = R, where L and R are the left-
 * and right-handed sides of the equations, given by 
 * L = 1 - P*V*G
 * R = P*V
 * Since G is expressed in an SWP basis, we also must include the basis-transormation matrices C */
void faddeev_dense_solver(cdouble*  U_array,
					      cdouble*  G_array,
					      int*		q_com_idx_array,	size_t num_q_com,
					      int*      deuteron_idx_array, size_t num_deuteron_states,
					      size_t    Nalpha,
					      size_t 	Nq_WP,
					      size_t 	Np_WP,
					      double**  CT_RM_array,
					      double**  VC_CM_array,
					      double*   P123_sparse_val_array,
					      int*      P123_sparse_row_array,
					      size_t*   P123_sparse_col_array,
					      size_t    P123_sparse_dim){
	
	/* Stores A and K arrays for the first indices in q_com_idx_array
	 * and deuteron_idx_array */
	bool store_first_A_array = false;
	bool store_first_K_array = false;

	/* Dense dimension of 3N-channel */
	size_t dense_dim = Nalpha * Nq_WP * Np_WP;
	
	std::complex<double>* L_array = new cdouble [dense_dim*dense_dim];
	std::complex<double>* R_array = new cdouble [dense_dim*dense_dim];

	double* CPVC_col_array 		   = new double [dense_dim];
	int*    CPVC_row_to_nnz_array  = new int    [dense_dim];
	int*    CPVC_nnz_to_row_array  = new int    [dense_dim];
	
	for (size_t j=0; j<num_q_com; j++){

		/* Reset L- and R-arrays */
		for (size_t idx=0; idx<dense_dim*dense_dim; idx++){
			L_array[idx] = 0;
			R_array[idx] = 0;
		}

		/* Construct L- and R-arrays */
		for (size_t idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
			for (size_t idx_q_c=0; idx_q_c<Nq_WP; idx_q_c++){
				for (size_t idx_p_c=0; idx_p_c<Np_WP; idx_p_c++){
					size_t col_idx = idx_alpha_c*Np_WP*Nq_WP + idx_q_c*Np_WP + idx_p_c;

					/* Reset CPVC-column array */
					for (size_t row_idx=0; row_idx<dense_dim; row_idx++){
						CPVC_col_array[row_idx] = 0;
						CPVC_row_to_nnz_array[row_idx] = -1;
						CPVC_nnz_to_row_array[row_idx] = -1;
					}

					/* Calculate CPVC-column */
					size_t CPVC_num_nnz = 0;
					calculate_CPVC_col(CPVC_col_array,
									   CPVC_row_to_nnz_array,
									   CPVC_nnz_to_row_array,
									   CPVC_num_nnz,
									   idx_alpha_c, idx_p_c, idx_q_c,
									   Nalpha, Nq_WP, Np_WP,
									   CT_RM_array,
									   VC_CM_array,
									   P123_sparse_val_array,
									   P123_sparse_row_array,
									   P123_sparse_col_array,
									   P123_sparse_dim);

    	    		//for (size_t row_idx=0; row_idx<dense_dim; row_idx++){
					int row_idx = 0;
					for (size_t nnz_idx=0; nnz_idx<CPVC_num_nnz; nnz_idx++){
						row_idx = CPVC_nnz_to_row_array[nnz_idx];

    	    		    //L_array[row_idx*dense_dim + col_idx] = -K_array[row*n + col];
						L_array[row_idx*dense_dim + col_idx] = -CPVC_col_array[nnz_idx]*G_array[j*dense_dim + col_idx];

    	    		    //R_array[row_idx*dense_dim + col_idx] =  A_array[row*n + col];
						R_array[row_idx*dense_dim + col_idx] =  CPVC_col_array[nnz_idx];
    	    		}

    	    		L_array[col_idx*dense_dim + col_idx] += 1;
				}
			}
    	}

		if (j==0 && store_first_A_array){
			store_array(R_array, dense_dim*dense_dim, "A_array.txt");
		}
		if (j==0 && store_first_K_array){
			store_array(L_array, dense_dim*dense_dim, "K_array.txt");
		}

		/* Solve */
		solve_MM(L_array, R_array, dense_dim);

		for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
			for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
				size_t idx_alpha_row = deuteron_idx_array[idx_d_row];
				size_t idx_alpha_col = deuteron_idx_array[idx_d_col];
				size_t idx_p_NDOS 	  = 0;
				size_t idx_q_NDOS 	  = q_com_idx_array[j];

				size_t idx_row_NDOS = idx_alpha_row*Nq_WP*Np_WP + idx_q_NDOS*Np_WP + idx_p_NDOS;
				size_t idx_col_NDOS = idx_alpha_col*Nq_WP*Np_WP + idx_q_NDOS*Np_WP + idx_p_NDOS;

				cdouble U_val = R_array[idx_row_NDOS*dense_dim + idx_col_NDOS];

				printf("   - U-matrix element for alpha'=%d, alpha=%d, q=%d: %.10e + %.10ei \n", idx_alpha_row, idx_alpha_col, idx_q_NDOS, idx_p_NDOS, U_val.real(), U_val.imag());
			}
		}

	}
	delete [] L_array;
	delete [] R_array;
	delete [] CPVC_col_array;
	delete [] CPVC_row_to_nnz_array;
	delete [] CPVC_nnz_to_row_array;
}

void pade_method_solve(cdouble*  U_array,
					   cdouble*  G_array,
					   int*		 q_com_idx_array,	 size_t num_q_com,
					   int*      deuteron_idx_array, size_t num_deuteron_states,
					   size_t    Nalpha,
					   size_t 	 Nq_WP,
					   size_t 	 Np_WP,
					   double**  CT_RM_array,
					   double**  VC_CM_array,
					   double*   P123_sparse_val_array,
					   int*      P123_sparse_row_array,
					   size_t*   P123_sparse_col_array,
					   size_t    P123_sparse_dim){

	/* Print Pade-approximant convergences */
	bool print_PA_convergences = false;
	/* Print Neumann terms */
	bool print_neumann_terms = true;
	
	/* Number of on-shell nucleon-deuteron channels (deuteron states can mix, hence ^2) */
	size_t num_on_shell_A_rows = num_deuteron_states * num_q_com;
	size_t num_on_shell_A_vals = num_deuteron_states * num_deuteron_states * num_q_com;
	
	/* Upper limit on polynomial approximation of Faddeev eq. */
	size_t N_pade = 14;
	size_t M_pade = 14;

	/* Coefficients for calculating Pade approximant */
	cdouble* a_coeff_array = new cdouble [ (N_pade + M_pade + 1) * num_on_shell_A_vals];

	/* Dense dimension of 3N-channel */
	size_t dense_dim = Nalpha * Nq_WP * Np_WP;

	/* Allocate row-arrays for A*A^n, where A=(C^T)(P)(VC) */
	cdouble* A_An_row_array 	 = new cdouble [dense_dim * num_on_shell_A_rows];
	cdouble* A_An_row_array_prev = new cdouble [dense_dim * num_on_shell_A_rows];

	/* Set A_An-arrays to zero */
	for (size_t i=0; i<dense_dim*num_on_shell_A_rows; i++){
		A_An_row_array[i] 	   = 0;
		A_An_row_array_prev[i] = 0;
	}
	
	//double* PVC_col_1 = new double [dense_dim];
	//double* PVC_col_2 = new double [dense_dim];
	//for (size_t i=0; i<dense_dim; i++){
	//	PVC_col_1[i] = 0;
	//	PVC_col_2[i] = 0;
	//}
	//calculate_PVC_col(PVC_col_1,
	//				  5, 0, 1,
	//				  Nalpha,      Nq_WP,   Np_WP,
	//				  VC_CM_array,
	//				  P123_sparse_val_array,
	//				  P123_sparse_row_array,
	//				  P123_sparse_col_array,
	//				  P123_sparse_dim);
	//calculate_PVC_col(PVC_col_2,
	//				  7, 0, 1,
	//				  Nalpha,      Nq_WP,   Np_WP,
	//				  VC_CM_array,
	//				  P123_sparse_val_array,
	//				  P123_sparse_row_array,
	//				  P123_sparse_col_array,
	//				  P123_sparse_dim);
	//for (size_t i=0; i<dense_dim; i++){
	//	if (PVC_col_1[i]!=0){
	//		std::cout << PVC_col_1[i] << std::endl;
	//		std::cout << PVC_col_2[i] << std::endl;
	//		std::cout << std::endl;
	//	}
	//}
	//raise_error("end of test");


	/* Set initial values for A_Kn_row_array, where K^n=1 for n=0 */
	printf("   - Starting CPVC-row calc. for all relevant rows \n");//row alpha=%d, q=%d, p=%d\n", idx_alpha_NDOS, idx_q_NDOS, idx_p_NDOS);
	auto timestamp_start = std::chrono::system_clock::now();
	/* Calculate CPVC-row and write to A_Kn_row_array_prev */
	calculate_all_CPVC_rows(A_An_row_array_prev,
							q_com_idx_array, num_q_com,
			   				deuteron_idx_array, num_deuteron_states,
							Nalpha, Nq_WP, Np_WP,
							CT_RM_array,
							VC_CM_array,
							P123_sparse_val_array,
							P123_sparse_row_array,
							P123_sparse_col_array,
							P123_sparse_dim);
	auto timestamp_end = std::chrono::system_clock::now();
	std::chrono::duration<double> time = timestamp_end - timestamp_start;
	printf("   - Time CPVC-rows: %.6f\n", time.count());

	/* First Neumann-term */
	printf("   - Pade iteration n=%d \n",0);
	for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
		for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
			for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
				/* Nucleon-deuteron on-shell (NDOS) indices
				 * (deuteron bound-state p-index is alwasy 0 due to eigenvalue ordering in SWP construction) */
				size_t idx_alpha_NDOS_row = deuteron_idx_array[idx_d_row];
				size_t idx_alpha_NDOS_col = deuteron_idx_array[idx_d_col];
				size_t idx_p_NDOS 	  = 0;
				size_t idx_q_NDOS 	  = q_com_idx_array[idx_q_com];

				size_t idx_row_NDOS   = idx_d_row*num_q_com + idx_q_com;
				size_t idx_col_NDOS   = idx_alpha_NDOS_col*Nq_WP*Np_WP + idx_q_NDOS*Np_WP + idx_p_NDOS;

				/* Calculate coefficient */
				cdouble a_coeff = A_An_row_array_prev[idx_row_NDOS*dense_dim + idx_col_NDOS];

				/* Store coefficient */
				size_t idx_NDOS = idx_d_row*num_deuteron_states*num_q_com + idx_d_col*num_q_com + idx_q_com;
				a_coeff_array[idx_NDOS*(N_pade+M_pade+1)] = a_coeff;
				
				if (print_neumann_terms){
					printf("     - Neumann term %d for alpha'=%d, alpha=%d, q=%d: %.16e + %.16ei \n", 0, idx_alpha_NDOS_row, idx_alpha_NDOS_col, idx_q_NDOS, a_coeff.real(), a_coeff.imag());
				}
			}
		}
	}
	
	/* Loop over number of Pade-terms we use */
	for (size_t n=1; n<N_pade+M_pade+1; n++){
		printf("   - Pade iteration n=%d \n",n);
		auto timestamp_start = std::chrono::system_clock::now();
		
		size_t counter_array [100];
		for (size_t i=0; i<100; i++){
			counter_array[i] = 0;
		}
		/* Iterate through columns of A */
		
		int num_threads = omp_get_max_threads();
		double*  times_array = new double [3*num_threads];
		for (int i=0; i<3*num_threads; i++){
			times_array[i] = 0;
		}

		#pragma omp parallel
		{
		/* Allocate row- and column-arrays for (C^T)(P)(VC) */
		double*  CPVC_col_array  		= new double [dense_dim];
		int*     CPVC_row_to_nnz_array  = new int    [dense_dim];
		int*     CPVC_nnz_to_row_array  = new int    [dense_dim];
		size_t   thread_idx = omp_get_thread_num();
		//cdouble* G_subarray		 = NULL;
		//cdouble* CPVCG_col_array = new cdouble [dense_dim];
		#pragma omp for
		for (size_t idx_q_c=0; idx_q_c<Nq_WP; idx_q_c++){
			for (size_t idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
				for (size_t idx_p_c=0; idx_p_c<Np_WP; idx_p_c++){
					size_t idx_col = idx_alpha_c*Nq_WP*Np_WP + idx_q_c*Np_WP + idx_p_c;

					double timestamp_0 = omp_get_wtime();
					/* Reset CPVC-column array */
					for (size_t row_idx=0; row_idx<dense_dim; row_idx++){
						CPVC_col_array[row_idx] = 0;

						CPVC_row_to_nnz_array[row_idx] = -1;
						CPVC_nnz_to_row_array[row_idx] = -1;
					}

					double timestamp_1 = omp_get_wtime();
					double time1 = timestamp_1 - timestamp_0;
					times_array[3*thread_idx] += time1;

					/* Calculate CPVC-column */
					size_t CPVC_num_nnz = 0;
					calculate_CPVC_col(CPVC_col_array,
									   CPVC_row_to_nnz_array,
									   CPVC_nnz_to_row_array,
									   CPVC_num_nnz,
									   idx_alpha_c, idx_p_c, idx_q_c,
									   Nalpha, Nq_WP, Np_WP,
									   CT_RM_array,
									   VC_CM_array,
									   P123_sparse_val_array,
									   P123_sparse_row_array,
									   P123_sparse_col_array,
									   P123_sparse_dim);
					counter_array[thread_idx] += CPVC_num_nnz;

					double timestamp_2 = omp_get_wtime();
					double time2 = timestamp_2 - timestamp_1;
					times_array[3*thread_idx + 1] += time2;
					

					/* Calculate all a-coefficients for calculated CPVC-column */
					for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
						for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
							size_t idx_row_NDOS = idx_d_row*num_q_com + idx_q_com;
							
							/* Dot product An*A */
							//cdouble inner_product = cdot_VV(&A_An_row_array_prev[idx_row_NDOS*dense_dim], CPVCG_col_array, dense_dim, 1, 1);

							cdouble inner_product = 0;
							////cdouble inner_products [omp_get_max_threads()];
							////#pragma omp parallel for
							//for (size_t i=0; i<dense_dim; i++){
							//	inner_product += A_An_row_array_prev[idx_row_NDOS*dense_dim + i] * CPVC_col_array[i] * G_array[idx_q_com*dense_dim + i];
							//	//inner_products[omp_get_thread_num()] += A_An_row_array_prev[idx_row_NDOS*dense_dim + i] * CPVC_col_array[i] * G_array[idx_q_com*dense_dim + i];
							//}

							int row_idx = 0;
							for (size_t nnz_idx=0; nnz_idx<CPVC_num_nnz; nnz_idx++){
								row_idx = CPVC_nnz_to_row_array[nnz_idx];
								inner_product += A_An_row_array_prev[idx_row_NDOS*dense_dim + row_idx] * CPVC_col_array[nnz_idx] * G_array[idx_q_com*dense_dim + row_idx];
							}

							A_An_row_array[idx_row_NDOS*dense_dim + idx_col] = inner_product;
						}
					}

					double timestamp_3 = omp_get_wtime();
					double time3 = timestamp_3 - timestamp_2;
					times_array[3*thread_idx + 2] += time3;
				}
			}
		}
		delete [] CPVC_col_array;
		delete [] CPVC_row_to_nnz_array;
		delete [] CPVC_nnz_to_row_array;
		}
		auto timestamp_end = std::chrono::system_clock::now();
		std::chrono::duration<double> time = timestamp_end - timestamp_start;

		double time_resetting_arrays = 0;
		double time_CPVC_cols = 0;
		double time_An_CPVC_multiply = 0;
		for (int thread_idx=0; thread_idx<num_threads; thread_idx++){
			time_resetting_arrays += times_array[3*thread_idx];
			time_CPVC_cols 		  += times_array[3*thread_idx+1];
			time_An_CPVC_multiply += times_array[3*thread_idx+2];
		}
		printf("     - Time resetting arrays:         %.6f \n", time_resetting_arrays/num_threads);
		printf("     - Time generating CPVC-cols:     %.6f \n", time_CPVC_cols/num_threads);
		printf("     - Time multiplying with An-rows: %.6f \n", time_An_CPVC_multiply/num_threads);
		printf("     - Total time:                    %.6f \n", time.count());
		printf("     - Done. \n");

		//size_t nnz_counts = 0;
		//for (size_t i=0; i<100; i++){
		//	nnz_counts += counter_array[i];
		//}
		//printf("NUMBER OF NNZ ELEMENTS IN A: %zu \n", nnz_counts);
		//printf("NUMBER OF NNZ ELEMENTS IN P: %zu \n", P123_sparse_dim);

		/* Rewrite previous A_An with current A_An */
		for (size_t i=0; i<num_on_shell_A_rows*dense_dim; i++){
			A_An_row_array_prev[i] = A_An_row_array[i];
		}

		printf("   - Extracting Neumann-series terms a_n \n");
		/* Extract coefficients "a" for Pade approximant */
		for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
			for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
				for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
					/* Nucleon-deuteron on-shell (NDOS) indices
					 * (deuteron bound-state p-index is alwasy 0 due to eigenvalue ordering in SWP construction) */
					size_t idx_alpha_NDOS_row = deuteron_idx_array[idx_d_row];
					size_t idx_alpha_NDOS_col = deuteron_idx_array[idx_d_col];
					size_t idx_p_NDOS 	  = 0;
					size_t idx_q_NDOS 	  = q_com_idx_array[idx_q_com];

					size_t idx_row_NDOS   = idx_d_row*num_q_com + idx_q_com;
					size_t idx_col_NDOS   = idx_alpha_NDOS_col*Nq_WP*Np_WP + idx_q_NDOS*Np_WP + idx_p_NDOS;
					
					/* Calculate coefficient */
					cdouble a_coeff = A_An_row_array[idx_row_NDOS*dense_dim + idx_col_NDOS];

					/* Store coefficient */
					size_t idx_NDOS = idx_d_row*num_deuteron_states*num_q_com + idx_d_col*num_q_com + idx_q_com;
					a_coeff_array[idx_NDOS*(N_pade+M_pade+1) + n] = a_coeff;

					if (print_neumann_terms){
						printf("     - Neumann term %d for alpha'=%d, alpha=%d, q=%d: %.16e + %.16ei \n", n, idx_alpha_NDOS_row, idx_alpha_NDOS_col, idx_q_NDOS, a_coeff.real(), a_coeff.imag());
					}
				}
			}
		}
		printf("     - Done \n");
	}

	/* Calculate Pade approximants (PA) */
	for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
		for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
			for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
				/* Nucleon-deuteron on-shell (NDOS) indices
				 * (deuteron bound-state p-index is alwasy 0 due to eigenvalue ordering in SWP construction) */
				size_t idx_alpha_NDOS_row = deuteron_idx_array[idx_d_row];
				size_t idx_alpha_NDOS_col = deuteron_idx_array[idx_d_col];
				size_t idx_q_NDOS 	   	  = q_com_idx_array[idx_q_com];

				size_t idx_NDOS = idx_d_row*num_deuteron_states*num_q_com + idx_d_col*num_q_com + idx_q_com;

				cdouble pade_approximants_array [N_pade+M_pade+1];

				double min_PA_diff = 1;
				size_t idx_best_PA = 0;

				/* Calculate PAs and find index for most converged PA */
				for (size_t NM=0; NM<N_pade+M_pade+1; NM+=2){
					cdouble PA = pade_approximant(&a_coeff_array[idx_NDOS*(N_pade+M_pade+1)], NM/2, NM/2, 1);

					pade_approximants_array[NM/2] = PA;

					double PA_diff = std::abs(PA - pade_approximants_array[NM/2 - 1]);
					if (print_PA_convergences){
						printf("PA[%d,%d] = %.16e + %.16ei, PA_diff = %.16e \n", NM/2,NM/2,PA.real(), PA.imag(), PA_diff);
					}
					/* Ignore PA_diff if numerically equal to the previous PA_diff, overwrite if smaller than min_PA_diff */
					if (PA_diff<min_PA_diff && PA_diff>1e-15){
						idx_best_PA = NM/2;
						min_PA_diff = PA_diff;
					}
				}

				/* Set U-matrix element equal "best" PA */
				U_array[idx_NDOS] = pade_approximants_array[idx_best_PA];
				printf(" - U-matrix element for alpha'=%d, alpha=%d, q=%d: %.10e + %.10ei \n", idx_alpha_NDOS_row, idx_alpha_NDOS_col, idx_q_NDOS, U_array[idx_NDOS].real(), U_array[idx_NDOS].imag());
			}
		}
	}

	delete [] A_An_row_array;
	delete [] A_An_row_array_prev;
	delete [] a_coeff_array;
}

void solve_faddeev_equations(cdouble*  U_array,
							 cdouble*  G_array,
							 double*   P123_sparse_val_array,
							 int*      P123_sparse_row_array,
							 int*      P123_sparse_col_array,
							 size_t    P123_sparse_dim,
							 double*   C_WP_unco_array,
							 double*   C_WP_coup_array,
							 double*   V_WP_unco_array,
							 double*   V_WP_coup_array,
							 int*	   q_com_idx_array,	   size_t num_q_com,
					   		 int*      deuteron_idx_array, size_t num_deuteron_states,
							 int       J_2N_max,
							 size_t    Nq_WP,
							 size_t    Np_WP,
							 size_t    Nalpha,
							 int*      L_2N_array,
							 int*      S_2N_array,
							 int*      J_2N_array,
							 int*      T_2N_array,
							 int*      L_1N_array, 
							 int*      two_J_1N_array){
	
	/* Test PVC- and CPVC-column multiplication routines with brute-force routines
	 * WARNING: VERY SLOW TEST, ONLY FOR BENCHMARKING */
	bool test_PVC_col_routine   = false;
	bool test_CPVC_col_routine  = false;

	/* Solve the Faddeev eq. using a dense MKL-solver.
	 * Obviously, this only works for small systems and is meant for benchmarking only */
	bool solve_dense = false;

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
									   T_2N_array,
									   L_1N_array,
									   two_J_1N_array);

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
									   T_2N_array,
									   L_1N_array,
									   two_J_1N_array);
	
	size_t  dense_dim = Nalpha * Nq_WP * Np_WP;

	//// Temp code to print P123 rows
	//size_t row_idx = 5*Np_WP*Nq_WP + 1*Np_WP + 0;
	//for (size_t nnz_idx=0; nnz_idx<P123_sparse_dim; nnz_idx++){
	//	if (P123_sparse_row_array[nnz_idx]==row_idx){
	//		std::cout << P123_sparse_val_array[nnz_idx] << std::endl;
	//	}
	//}
	
	/* Convert row-major sparse format to column-major */
	printf(" - Converting P123 from row- to column-major ... \n");
	unsorted_sparse_to_coo_col_major_sorter(&P123_sparse_val_array,
											&P123_sparse_row_array,
											&P123_sparse_col_array,
											P123_sparse_dim,
											dense_dim);
	printf("   - Done \n");

	//// Temp code to print P123 cols (easier after RM to CM change)
	//size_t col_idx = 5*Np_WP*Nq_WP + 1*Np_WP + 0;
	//for (size_t nnz_idx=0; nnz_idx<P123_sparse_dim; nnz_idx++){
	//	if (P123_sparse_col_array[nnz_idx]==col_idx){
	//		std::cout << P123_sparse_val_array[nnz_idx] << std::endl;
	//	}
	//}

	/* Convert from COO format to CSC format */
	printf(" - Converting P123 from COO to CSC format ... \n");
	size_t* P123_sparse_col_array_csc = new size_t [dense_dim];
	coo_to_csr_format_converter(P123_sparse_col_array,
								P123_sparse_col_array_csc,
								P123_sparse_dim,
								dense_dim);
	printf("   - Done \n");
	
	/* Test optimized routine for PVC columns */
	if (test_PVC_col_routine){
		printf(" - Testing PVC-column routine ... \n");
		PVC_col_calc_test(Nalpha,
						  Nq_WP,
						  Np_WP,
						  VC_CM_array,
						  P123_sparse_val_array,
						  P123_sparse_row_array,
						  P123_sparse_col_array_csc,
						  P123_sparse_dim);
		printf("   - Done \n");
	}
	/* Test optimized routine for CPVC columns */
	if (test_CPVC_col_routine){
		printf(" - Testing CPVC-column routine ... \n");
		CPVC_col_calc_test(Nalpha,
						   Nq_WP,
						   Np_WP,
						   CT_RM_array,
						   VC_CM_array,
						   P123_sparse_val_array,
						   P123_sparse_row_array,
						   P123_sparse_col_array_csc,
						   P123_sparse_dim);
		printf("   - Done \n");
	}
	
	printf(" - Solving Faddeev equation ... \n");
	auto timestamp_solve_start = std::chrono::system_clock::now();
	
	pade_method_solve(U_array,
					  G_array,
					  q_com_idx_array,	  num_q_com,
					  deuteron_idx_array, num_deuteron_states,
					  Nalpha,
					  Nq_WP,
					  Np_WP,
					  CT_RM_array,
					  VC_CM_array,
					  P123_sparse_val_array,
					  P123_sparse_row_array,
					  P123_sparse_col_array_csc,
					  P123_sparse_dim);

	auto timestamp_solve_end = std::chrono::system_clock::now();
	std::chrono::duration<double> time_solve = timestamp_solve_end - timestamp_solve_start;
	printf("   - Done. Time used: %.6f\n", time_solve.count());

	std::string U_mat_foldername = "../../Data/U_matrix_elements/";
	std::string U_mat_filename = U_mat_foldername + "U_PW_elements_Np_" + std::to_string(Np_WP)
																		+ "_Nq_" + std::to_string(Np_WP)
																		+ "_JP_" + std::to_string(3)
																		+ "_" + std::to_string(1)
																		+ "_Jmax_" + std::to_string(1)
																		+ ".csv";
	store_U_matrix_elements_csv(U_array,
							    q_com_idx_array,    (size_t) num_q_com,
					  		    deuteron_idx_array, (size_t) num_deuteron_states,
							    L_1N_array, 
							    two_J_1N_array,
							    U_mat_filename);
	
	if (solve_dense){
		printf(" - Solving Faddeev equation using a dense direct solver (WARNING: CAN TAKE LONG) ... \n");
		auto timestamp_solve_start = std::chrono::system_clock::now();

		faddeev_dense_solver(U_array,
						     G_array,
						     q_com_idx_array,	 num_q_com,
						     deuteron_idx_array, num_deuteron_states,
						     Nalpha,
						     Nq_WP,
						     Np_WP,
						     CT_RM_array,
						     VC_CM_array,
						     P123_sparse_val_array,
						     P123_sparse_row_array,
						     P123_sparse_col_array_csc,
						     P123_sparse_dim);

		auto timestamp_solve_end = std::chrono::system_clock::now();
		std::chrono::duration<double> time_solve = timestamp_solve_end - timestamp_solve_start;
		printf("   - Done. Time used: %.6f\n", time_solve.count());
	}
}

/* Create array of pointers to C^T matrices for product (C^T)PVC in row-major format */
void create_CT_row_maj_3N_pointer_array(double** CT_RM_array,
										double*  C_WP_unco_array,
										double*  C_WP_coup_array,
										size_t   Np_WP,
										int      J_2N_max,
										size_t   Nalpha,
										int*     L_2N_array,
										int*     S_2N_array,
										int*     J_2N_array,
										int*     T_2N_array,
							 			int*     L_1N_array, 
							 			int*     two_J_1N_array){

	double* C_subarray  = NULL;
	double* CT_subarray = NULL;

	double* CT_unco_array = new double [Np_WP*Np_WP   * 2*(J_2N_max+1)];
	double* CT_coup_array = new double [Np_WP*Np_WP*4 *    J_2N_max   ];
	
	/* Copy and transpose all 2N-uncoupled C-arrays */
	for (int J_2N=0; J_2N<J_2N_max+1; J_2N++){
		for (int S_2N=0; S_2N<2; S_2N++){
			size_t idx_chn_unco       = 2*J_2N + S_2N;
			size_t idx_2N_mat_WP_unco = idx_chn_unco*Np_WP*Np_WP;
			
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
		size_t idx_chn_coup     = J_2N-1;
		size_t idx_2N_mat_WP_coup = idx_chn_coup*4*Np_WP*Np_WP;

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
	for (size_t idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
		int L_2N_r 	   = L_2N_array[idx_alpha_r];
		int S_2N_r 	   = S_2N_array[idx_alpha_r];
		int J_2N_r 	   = J_2N_array[idx_alpha_r];
		int T_2N_r 	   = T_2N_array[idx_alpha_r];
		int L_1N_r 	   = L_1N_array[idx_alpha_r];
		int two_J_1N_r = two_J_1N_array[idx_alpha_r];

		/* Column state */
		for (size_t idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
			int L_2N_c 	   = L_2N_array[idx_alpha_c];
			int S_2N_c 	   = S_2N_array[idx_alpha_c];
			int J_2N_c 	   = J_2N_array[idx_alpha_c];
			int T_2N_c 	   = T_2N_array[idx_alpha_c];
			int L_1N_c 	   = L_1N_array[idx_alpha_c];
			int two_J_1N_c = two_J_1N_array[idx_alpha_c];

			/* Check if possible channel through interaction */
			if (T_2N_r==T_2N_c &&
			    J_2N_r==J_2N_c &&
				S_2N_r==S_2N_c &&
				abs(L_2N_r-L_2N_c)<=2 &&
				L_1N_r==L_1N_c &&
				two_J_1N_r==two_J_1N_c){

				/* Detemine if this is a coupled channel */
				bool coupled_matrix = false;
				if (L_2N_r!=L_2N_c || (L_2N_r==L_2N_c & L_2N_r!=J_2N_r & J_2N_r!=0)){ // This counts 3P0 as uncoupled; used in matrix structure
					coupled_matrix  = true;
				}

				/* find which VC-product corresponds to the current coupling */
				if (coupled_matrix){
					size_t idx_chn_coup       = (size_t) J_2N_r-1;
					size_t idx_2N_mat_WP_coup = idx_chn_coup*4*Np_WP*Np_WP;
					if (L_2N_r<L_2N_c){       // L_r=J_r-1, L_c=J_r+1
						CT_subarray = &CT_coup_array[idx_2N_mat_WP_coup + 1*Np_WP*Np_WP];
					}
					else if (L_2N_r>L_2N_c){  // L_r=J_r+1, L_c=J_r-1
						CT_subarray = &CT_coup_array[idx_2N_mat_WP_coup + 2*Np_WP*Np_WP];
					}
					else if (L_2N_r<J_2N_c){  // L_r=J_r-1, L_c=J_r-1
						CT_subarray = &CT_coup_array[idx_2N_mat_WP_coup + 0*Np_WP*Np_WP];
					}
					else{               // L_r=J_r+1, L_c=J_r+1
						CT_subarray = &CT_coup_array[idx_2N_mat_WP_coup + 3*Np_WP*Np_WP];
					}
				}
				else{
					size_t idx_chn_unco       = (size_t) 2*J_2N_r + S_2N_r;
					size_t idx_2N_mat_WP_unco = idx_chn_unco*Np_WP*Np_WP;
					CT_subarray = &CT_unco_array[idx_2N_mat_WP_unco];
				}
			}
			else{
				CT_subarray = NULL;
			}

			/* Unique index for two given states alpha */
			size_t idx_CT_RM = idx_alpha_r*Nalpha + idx_alpha_c;
			CT_RM_array[idx_CT_RM] = CT_subarray;
		}
	}
	
}

/* Restructures NN coupled matrix as 4 seperate matrices
 * !!! WARNING: COLUMN-MAJOR ALGORITHM !!! */
void restructure_coupled_VC_product(double* VC_product, size_t Np_WP){
	double* VC_product_temp = new double [4*Np_WP*Np_WP];

	for (size_t col_block=0; col_block<2; col_block++){
		for (size_t row_block=0; row_block<2; row_block++){

			size_t idx_block = col_block*2*Np_WP*Np_WP + row_block*Np_WP*Np_WP;

			for (size_t col=0; col<Np_WP; col++){
				for (size_t row=0; row<Np_WP; row++){
					size_t prestructure_col_idx = (col + Np_WP*col_block)*2*Np_WP;
					size_t prestructure_row_idx =  row + Np_WP*row_block;

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
										size_t   Np_WP,
										int      J_2N_max,
										size_t   Nalpha,
										int*     L_2N_array,
										int*     S_2N_array,
										int*     J_2N_array,
										int*     T_2N_array,
							 			int*     L_1N_array, 
							 			int*     two_J_1N_array){

	double* V_subarray = NULL;
	double* C_subarray = NULL;

	double* VC_unco_array = new double [Np_WP*Np_WP   * 2*(J_2N_max+1)];
	double* VC_coup_array = new double [Np_WP*Np_WP*4 *    J_2N_max   ];

	double* VC_product = NULL;

	/* Calculate all 2N-uncoupled VC-products and convert to column-major format */
	for (int J_2N=0; J_2N<J_2N_max+1; J_2N++){
		for (int S_2N=0; S_2N<2; S_2N++){
			size_t idx_chn_unco       = (size_t) 2*J_2N + S_2N;
			size_t idx_2N_mat_WP_unco = idx_chn_unco*Np_WP*Np_WP;
			
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
		size_t idx_chn_coup       = (size_t) J_2N-1;
		size_t idx_2N_mat_WP_coup = idx_chn_coup*4*Np_WP*Np_WP;

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
	for (size_t idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
		int L_2N_r     = L_2N_array[idx_alpha_r];
		int S_2N_r     = S_2N_array[idx_alpha_r];
		int J_2N_r     = J_2N_array[idx_alpha_r];
		int T_2N_r     = T_2N_array[idx_alpha_r];
		int L_1N_r 	   = L_1N_array[idx_alpha_r];
		int two_J_1N_r = two_J_1N_array[idx_alpha_r];

		/* Column state */
		for (size_t idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
			int L_2N_c     = L_2N_array[idx_alpha_c];
			int S_2N_c     = S_2N_array[idx_alpha_c];
			int J_2N_c     = J_2N_array[idx_alpha_c];
			int T_2N_c     = T_2N_array[idx_alpha_c];
			int L_1N_c 	   = L_1N_array[idx_alpha_c];
			int two_J_1N_c = two_J_1N_array[idx_alpha_c];

			/* Check if possible channel through interaction */
			if (T_2N_r==T_2N_c &&
			    J_2N_r==J_2N_c &&
				S_2N_r==S_2N_c &&
				abs(L_2N_r-L_2N_c)<=2 &&
				L_1N_r==L_1N_c &&
				two_J_1N_r==two_J_1N_c){

				/* Detemine if this is a coupled channel */
				bool coupled_matrix = false;
				if (L_2N_r!=L_2N_c || (L_2N_r==L_2N_c && L_2N_r!=J_2N_r && J_2N_r!=0)){ // This counts 3P0 as uncoupled; used in matrix structure
					coupled_matrix  = true;
				}

				/* find which VC-product corresponds to the current coupling */
				if (coupled_matrix){
					size_t idx_chn_coup       = (size_t) J_2N_r-1;
					size_t idx_2N_mat_WP_coup = idx_chn_coup*4*Np_WP*Np_WP;
					if (L_2N_r<L_2N_c){       // L_r=J_r-1, L_c=J_r+1
						VC_product = &VC_coup_array[idx_2N_mat_WP_coup + 1*Np_WP*Np_WP];
					}
					else if (L_2N_r>L_2N_c){  // L_r=J_r+1, L_c=J_r-1
						VC_product = &VC_coup_array[idx_2N_mat_WP_coup + 2*Np_WP*Np_WP];
					}
					else if (L_2N_r<J_2N_c){  // L_r=J_r-1, L_c=J_r-1
						VC_product = &VC_coup_array[idx_2N_mat_WP_coup + 0*Np_WP*Np_WP];
					}
					else{               // L_r=J_r+1, L_c=J_r+1
						VC_product = &VC_coup_array[idx_2N_mat_WP_coup + 3*Np_WP*Np_WP];
					}
				}
				else{
					size_t idx_chn_unco       = (size_t) 2*J_2N_r + S_2N_r;
					size_t idx_2N_mat_WP_unco = idx_chn_unco*Np_WP*Np_WP;
					VC_product = &VC_unco_array[idx_2N_mat_WP_unco];
				}
			}
			else{
				VC_product = NULL;
			}

			/* Unique index for two given states alpha */
			size_t idx_VC_CM = idx_alpha_r*Nalpha + idx_alpha_c;
			VC_CM_array[idx_VC_CM] = VC_product;
		}
	}
}

void PVC_col_brute_force(double*  col_array,
					     size_t   idx_alpha_c, size_t idx_p_c, size_t idx_q_c,
					     size_t   Nalpha,      size_t Nq_WP,   size_t Np_WP,
					     double** VC_CM_array,
					     double*  P123_val_array,
					     int*  	  P123_row_array,
					     size_t*  P123_col_array,
					     size_t   P123_dim){
	
	bool   print_content = false;
	size_t idx1 = idx_alpha_c*Np_WP*Nq_WP + idx_q_c*Np_WP + idx_p_c;
	double* VC_ptr     = NULL;

	for (size_t idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
		for (size_t idx_q_r=0; idx_q_r<Nq_WP; idx_q_r++){
			for (size_t idx_p_r=0; idx_p_r<Np_WP; idx_p_r++){
				size_t idx_P123_row = idx_alpha_r*Nq_WP*Np_WP + idx_q_r*Np_WP + idx_p_r;

				double sum = 0; 
				for (size_t idx_alpha_j=0; idx_alpha_j<Nalpha; idx_alpha_j++){
					VC_ptr = VC_CM_array[idx_alpha_c*Nalpha + idx_alpha_j];
					if (VC_ptr!=NULL){
						for (size_t idx_q_j=0; idx_q_j<Nq_WP; idx_q_j++){
							if (idx_q_j == idx_q_c){
								for (size_t idx_p_j=0; idx_p_j<Np_WP; idx_p_j++){
									if (print_content){std::cout << "  - Index j " << idx_alpha_j << " " << idx_q_j << " " << idx_p_j << std::endl;}

									size_t idx_P123_col = idx_alpha_j*Nq_WP*Np_WP + idx_q_j*Np_WP + idx_p_j;

									size_t idx_P123_row_lower = P123_col_array[idx_P123_col];
									size_t idx_P123_row_upper = P123_col_array[idx_P123_col+1];
									bool idx_found = false;
									size_t idx_P123_val = 0;
									for (size_t nnz_idx=idx_P123_row_lower; nnz_idx<idx_P123_row_upper; nnz_idx++){
										if (P123_row_array[nnz_idx] == idx_P123_row){
											idx_found = true;
											idx_P123_val = nnz_idx;
										}
									}

									if (idx_found){
										double P_element  = P123_val_array[idx_P123_val];
										double VC_element = VC_ptr[idx_p_c*Np_WP + idx_p_j];

										sum += P_element * VC_element;
									}
								}
							}
						}
					}
				}

				size_t idx_row = idx_alpha_r*Np_WP*Nq_WP + idx_q_r*Np_WP + idx_p_r;
				/* Write inner product to col_array of PVC-product */
				col_array[idx_row] = sum;
			}
		}
	}
}

void CPVC_col_brute_force(double*  col_array,
						  size_t   idx_alpha_c, size_t idx_p_c, size_t idx_q_c,
						  size_t   Nalpha,      size_t Nq_WP,   size_t Np_WP,
						  double** CT_RM_array,
						  double** VC_CM_array,
						  double*  P123_val_array,
						  int*     P123_row_array,
						  size_t*  P123_col_array,
						  size_t   P123_dim){
	double* CT_ptr = NULL;
	double* VC_ptr = NULL;

	bool print_content = false;

	/* Generate PVC-column */
	double* PVC_col = new double [Nalpha*Nq_WP*Np_WP];
	/* Ensure PVC_col contains only zeroes */
	for (size_t idx=0; idx<Nalpha*Nq_WP*Np_WP; idx++){
		PVC_col[idx] = 0;
	}
	
	PVC_col_brute_force(PVC_col,
					    idx_alpha_c, idx_p_c, idx_q_c,
					    Nalpha,      Nq_WP,   Np_WP,
					    VC_CM_array,
					    P123_val_array,
					    P123_row_array,
					    P123_col_array,
					    P123_dim);

	/* Index: r */
	for (size_t idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
		for (size_t idx_q_r=0; idx_q_r<Nq_WP; idx_q_r++){
			for (size_t idx_p_r=0; idx_p_r<Np_WP; idx_p_r++){
				if (print_content){std::cout << "Index r " << idx_alpha_r << " " << idx_q_r << " " << idx_p_r << std::endl;}
				double sum = 0;
				/* Index: i */
				for (size_t idx_alpha_i=0; idx_alpha_i<Nalpha; idx_alpha_i++){
					CT_ptr = CT_RM_array[idx_alpha_r*Nalpha + idx_alpha_i];
					if (CT_ptr!=NULL){
						for (size_t idx_q_i=0; idx_q_i<Nq_WP; idx_q_i++){
							if (idx_q_r == idx_q_i){
								for (size_t idx_p_i=0; idx_p_i<Np_WP; idx_p_i++){
									if (print_content){std::cout << " - Index i " << idx_alpha_i << " " << idx_q_i << " " << idx_p_i << std::endl;}

									size_t PVC_col_idx = idx_alpha_i*Np_WP*Nq_WP + idx_q_i*Np_WP + idx_p_i;

									double CT_element  = CT_ptr[idx_p_r*Np_WP + idx_p_i];
									double PVC_element = PVC_col[PVC_col_idx];
									sum += CT_element * PVC_element;
								}
							}
						}
					}
				}

				size_t idx_CPVC = idx_alpha_r*Nq_WP*Np_WP + idx_q_r*Np_WP + idx_p_r;
				col_array[idx_CPVC] = sum;
			}
		}
	}
}

void PVC_col_calc_test(size_t   Nalpha,
					   size_t 	Nq_WP,
					   size_t 	Np_WP,
					   double** VC_CM_array,
					   double*  P123_sparse_val_array,
					   int*  	P123_sparse_row_array,
					   size_t*  P123_sparse_col_array,
					   size_t   P123_sparse_dim){

	bool print_content = false;

	size_t dense_dim = Nalpha * Nq_WP * Np_WP;

	/* Create column of (C^T)(P)(VC) */
	double* PVC_col_array    = new double [dense_dim];
	double* PVC_col_array_BF = new double [dense_dim];

	bool all_cols_are_zero = true;
	for (size_t idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
		for (size_t idx_q_c=0; idx_q_c<Nq_WP; idx_q_c++){
			for (size_t idx_p_c=0; idx_p_c<Np_WP; idx_p_c++){
				
				if (idx_alpha_c<=2){continue;}

				/* Reset array */
				for (size_t idx=0; idx<dense_dim; idx++){
					PVC_col_array[idx]    = 0;
					PVC_col_array_BF[idx] = 0;
				}
	
				/* Calculate CPVC-column */
				calculate_PVC_col(PVC_col_array,
								  idx_alpha_c, idx_p_c, idx_q_c,
								  Nalpha, Nq_WP, Np_WP,
								  VC_CM_array,
								  P123_sparse_val_array,
								  P123_sparse_row_array,
								  P123_sparse_col_array,
								  P123_sparse_dim);
				
				PVC_col_brute_force(PVC_col_array_BF,
									idx_alpha_c, idx_p_c, idx_q_c,
									Nalpha, Nq_WP, Np_WP,
									VC_CM_array,
									P123_sparse_val_array,
									P123_sparse_row_array,
									P123_sparse_col_array,
									P123_sparse_dim);
				
				double col_sum =0;
				for (size_t idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
					for (size_t idx_q_r=0; idx_q_r<Nq_WP; idx_q_r++){
						for (size_t idx_p_r=0; idx_p_r<Np_WP; idx_p_r++){
							
							size_t idx_c = idx_alpha_c*Np_WP*Nq_WP + idx_q_c*Np_WP + idx_p_c;
							size_t idx_r = idx_alpha_r*Np_WP*Nq_WP + idx_q_r*Np_WP + idx_p_r;
							double diff = abs(PVC_col_array[idx_r] - PVC_col_array_BF[idx_r]);

							if ( diff > 1e-15 ){
								printf("   - Discrepency found: \n");
								printf("     - Row %d: alpha=%d, q=%d, p=%d \n", idx_r, idx_alpha_r, idx_q_r, idx_p_r);
								printf("     - Col %d: alpha=%d, q=%d, p=%d \n", idx_c, idx_alpha_c, idx_q_c, idx_p_c);
								printf("     - Discrepency: %.16f (optimized) vs. %.16f (brute force)\n", PVC_col_array[idx_r], PVC_col_array_BF[idx_r]);
								raise_error("PVC benchmarking failed");
							}
							col_sum += PVC_col_array[idx_r];
						}
					}
				}
				if (col_sum!=0 ){
					all_cols_are_zero = false;
				}
			}
		}
	}
	if (all_cols_are_zero){
		printf("   - All PVC-columns are zero! PVC-test likely failed \n");
	}
	
	delete [] PVC_col_array;
	delete [] PVC_col_array_BF;	
}

void CPVC_col_calc_test(size_t   Nalpha,
						size_t 	 Nq_WP,
						size_t 	 Np_WP,
						double** CT_RM_array,
						double** VC_CM_array,
						double*  P123_sparse_val_array,
						int*     P123_sparse_row_array,
						size_t*  P123_sparse_col_array,
						size_t   P123_sparse_dim){

	bool print_content = false;

	size_t dense_dim = Nalpha * Nq_WP * Np_WP;

	/* Create column of (C^T)(P)(VC) */
	double* CPVC_col_array    	   = new double [dense_dim];
	int*    CPVC_row_to_nnz_array  = new int    [dense_dim];
	int*    CPVC_nnz_to_row_array  = new int    [dense_dim];
	double* CPVC_col_array_BF 	   = new double [dense_dim];

	bool all_cols_are_zero = true;
	for (size_t idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
		for (size_t idx_q_c=0; idx_q_c<Nq_WP; idx_q_c++){
			for (size_t idx_p_c=0; idx_p_c<Np_WP; idx_p_c++){
	
				/* Reset array */
				for (size_t idx=0; idx<dense_dim; idx++){
					CPVC_col_array[idx]    	   = 0;
					CPVC_row_to_nnz_array[idx] = -1;
					CPVC_nnz_to_row_array[idx] = -1;
					CPVC_col_array_BF[idx]	   = 0;
				}
	
				/* Calculate CPVC-column */
				size_t CPVC_num_nnz = 0;
				calculate_CPVC_col(CPVC_col_array,
								   CPVC_row_to_nnz_array,
								   CPVC_nnz_to_row_array,
								   CPVC_num_nnz,
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
				
				double col_sum =0;
				for (size_t idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
					for (size_t idx_q_r=0; idx_q_r<Nq_WP; idx_q_r++){
						for (size_t idx_p_r=0; idx_p_r<Np_WP; idx_p_r++){
							
							size_t idx_c = idx_alpha_c*Np_WP*Nq_WP + idx_q_c*Np_WP + idx_p_c;
							size_t idx_r = idx_alpha_r*Np_WP*Nq_WP + idx_q_r*Np_WP + idx_p_r;

							size_t nnz_idx = CPVC_row_to_nnz_array[idx_r];

							double diff = abs(CPVC_col_array[nnz_idx] - CPVC_col_array_BF[idx_r]);

							if ( diff > 1e-15 ){
								printf("   - Discrepency found: \n");
								printf("     - Row %d: alpha=%d, q=%d, p=%d \n", idx_r, idx_alpha_r, idx_q_r, idx_p_r);
								printf("     - Col %d: alpha=%d, q=%d, p=%d \n", idx_c, idx_alpha_c, idx_q_c, idx_p_c);
								printf("     - Discrepency: %.16f (optimized) vs. %.16f (brute force)\n", CPVC_col_array[nnz_idx], CPVC_col_array_BF[idx_r]);
								raise_error("CPVC benchmarking failed");
							}
							col_sum += CPVC_col_array[idx_r];
						}
					}
				}
				if (col_sum!=0 ){
					all_cols_are_zero = false;
				}
			}
		}
	}
	if (all_cols_are_zero){
		printf("   - All CPVC-columns are zero! PVC-test likely failed \n");
	}
	
	delete [] CPVC_col_array;
	delete [] CPVC_row_to_nnz_array;
	delete [] CPVC_nnz_to_row_array;
	delete [] CPVC_col_array_BF;
}

/* Solves Faddeev on the form
 * (1-AG)U = A
 * where A = C^T PVC */
//void direct_sparse_solve(cdouble*  U_array,
//						 cdouble*  G_array,
//						 int       idx_on_shell,
//						 int       Nalpha,
//						 int 	   Nq_WP,
//						 int 	   Np_WP,
//						 double**  CT_RM_array,
//						 double**  VC_CM_array,
//						 double*   P123_sparse_val_array,
//						 int*      P123_sparse_row_array,
//						 int*      P123_sparse_col_array,
//						 int       P123_sparse_dim){
//	
//	/* Dense dimension of 3N-channel */
//	size_t dense_dim = Nalpha * Nq_WP * Np_WP;
//
//	/* Sparse dimension and step-length for incrementing sparse array size */
//	int		A_sparse_dim	   = 0;
//	int 	sparse_step_length = 0;
//	if (dense_dim>1000){
//		sparse_step_length = dense_dim/1000;
//	}
//	else{
//		sparse_step_length = dense_dim;
//	}
//	int current_array_dim  = sparse_step_length;
//
//	/* Dynamically sized sparse-storage COO (CM) array format for A-matrix */
//	double* A_sparse_val_array = new double [sparse_step_length];
//	int* 	A_sparse_row_array = new int    [sparse_step_length];
//	int* 	A_sparse_col_array = new int    [sparse_step_length];
//
//	/* Allocate column-array for (C^T)(P)(VC) */
//	double 				 CPVC_col_array 	  [dense_dim];
//
//	std::complex<double> A_on_shell_col_array [dense_dim];
//	std::complex<double>* A_dense_array = new std::complex<double> [dense_dim*dense_dim];
//	//double A_on_shell_col_array [dense_dim];
//	//double* A_dense_array	   = new double [dense_dim*dense_dim];
//	
//	for (int idx=0; idx<dense_dim*dense_dim; idx++){
//		A_dense_array[idx] = 0;
//	}
//
//	/* Calculation of columns of A */
//	for (int idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
//		for (int idx_q_c=0; idx_q_c<Nq_WP; idx_q_c++){
//			for (int idx_p_c=0; idx_p_c<Np_WP; idx_p_c++){
//
//				/* Reset CPVC-column array */
//				for (int row_idx=0; row_idx<dense_dim; row_idx++){
//					CPVC_col_array[row_idx]    = 0;
//				}
//
//				/* Calculate CPVC-column */
//				calculate_CPVC_col(CPVC_col_array,
//								   idx_alpha_c, idx_p_c, idx_q_c,
//								   Nalpha, Nq_WP, Np_WP,
//								   CT_RM_array,
//								   VC_CM_array,
//								   P123_sparse_val_array,
//								   P123_sparse_row_array,
//								   P123_sparse_col_array,
//								   P123_sparse_dim);
//				
//				/* Append on-shell column to right-hand side vector of Faddeev eq. (and convert type to complex) */
//				int col_idx = idx_alpha_c*Np_WP*Nq_WP + idx_q_c*Np_WP + idx_p_c;
//				if (col_idx == idx_on_shell){
//					for (int row_idx=0; row_idx<dense_dim; row_idx++){
//						A_on_shell_col_array[row_idx] = CPVC_col_array[row_idx];
//					}
//				}
//				
//				/* Loop through rows of CPVC-column and append to A if non-zero */
//				for (int row_idx=0; row_idx<dense_dim; row_idx++){
//
//					double element = CPVC_col_array[row_idx];
//					/* Add element if non-zero or if this a diagonal element (needed for "1-AG" in Faddeev eq.) */
//					if (element!=0 || col_idx==row_idx){
//
//						/* Append to sparse value and index arrays */
//						A_sparse_val_array[A_sparse_dim] = element;
//						A_sparse_row_array[A_sparse_dim] = row_idx;
//						A_sparse_col_array[A_sparse_dim] = col_idx;
//
//						A_dense_array[row_idx*dense_dim + col_idx] = element;
//
//						/* Increment sparse dimension (num of non-zero elements) */
//						A_sparse_dim += 1;
//
//						/* If the dimension goes over the array dimension we increase the array size
//						 * via a copy-paste-type routine, and increment the current array dimension */
//						if ( A_sparse_dim>=current_array_dim ){  // This should occur a small amount of the time
//							increase_sparse_array_size(&A_sparse_val_array, current_array_dim, sparse_step_length);
//							increase_sparse_array_size(&A_sparse_row_array, current_array_dim, sparse_step_length);
//							increase_sparse_array_size(&A_sparse_col_array, current_array_dim, sparse_step_length);
//
//							/* Increment sparse-array dimension */
//							current_array_dim += sparse_step_length;
//						}
//					}
//				}
//			}
//		}
//	}
//	/* Contract arrays to minimal size (number of non-zero elements) */
//	reduce_sparse_array_size(&A_sparse_val_array, current_array_dim, A_sparse_dim);
//	reduce_sparse_array_size(&A_sparse_row_array, current_array_dim, A_sparse_dim);
//	reduce_sparse_array_size(&A_sparse_col_array, current_array_dim, A_sparse_dim);
//
//	/* Conversion from column-major COO to row-major COO array format for A-matrix 
//	 * !!! WARNING: THIS ROUTINE REWRITES INPUT ARRAYS TO NEW FORMAT, OLD FORMAT IS DELETED !!! */
//	coo_col_major_to_coo_row_major_converter(&A_sparse_val_array,
//											 &A_sparse_row_array,
//											 &A_sparse_col_array,
//											 A_sparse_dim,
//											 dense_dim);
//	/* Conversion from COO to CSR array format for A-matrix */
//	int* A_idx_row_array_csr = new int [dense_dim+1];
//	coo_to_csr_format_converter(A_sparse_row_array,
//                                A_idx_row_array_csr,
//                                A_sparse_dim,
//                                dense_dim);
//
//	/* Test sparse format */
//	if (false){
//		for (int row_idx=0; row_idx<dense_dim; row_idx++){
//			int nnz_idx_lower = A_idx_row_array_csr[row_idx];
//			int nnz_idx_upper = A_idx_row_array_csr[row_idx+1];
//			for (int nnz_idx=nnz_idx_lower; nnz_idx<nnz_idx_upper; nnz_idx++){
//				int col_idx = A_sparse_col_array[nnz_idx];
//
//				//double val_sparse = A_sparse_val_array[nnz_idx];
//				//double val_dense  = A_dense_array[row_idx*dense_dim + col_idx];
//				std::complex<double> val_sparse = A_sparse_val_array[nnz_idx];
//				std::complex<double> val_dense  = A_dense_array[row_idx*dense_dim + col_idx];
//				if (val_sparse!=val_dense){
//					std::cout << "Mismatch. Sparse val: " << val_sparse.real() << " " << val_sparse.imag() << std::endl;
//					std::cout << "Mismatch. Dense val:  " << val_dense.real() << " " << val_dense.imag() << std::endl;
//					raise_error("Sparse format conversion failed");
//				}
//			}
//		}
//	}
//
//	/* Identity "I" minus "AG"-product and conversion to complex type */
//	std::complex<double>* IAG_sparse_val_array_cmplx = new std::complex<double> [A_sparse_dim];
//	std::complex<double>* IAG_dense_array_cmplx = new std::complex<double> [dense_dim*dense_dim];
//	//double* IAG_sparse_val_array_cmplx = new double [A_sparse_dim];
//	//double* IAG_dense_array_cmplx	   = new double [dense_dim*dense_dim];
//	for (int idx=0; idx<dense_dim*dense_dim; idx++){
//		IAG_dense_array_cmplx[idx] = 0;
//	}
//	for (int row_idx=0; row_idx<dense_dim; row_idx++){
//
//		int nnz_idx_lower = A_idx_row_array_csr[row_idx];
//		int nnz_idx_upper = A_idx_row_array_csr[row_idx+1];
//
//		for (int nnz_idx=nnz_idx_lower; nnz_idx<nnz_idx_upper; nnz_idx++){
//
//			int col_idx = A_sparse_col_array[nnz_idx];
//
//			std::complex<double> G_val = G_array[col_idx];
//
//			/* AG-multiplication and minus-sign */
//			IAG_sparse_val_array_cmplx[nnz_idx] = -A_sparse_val_array[nnz_idx] * G_val;
//			//IAG_sparse_val_array_cmplx[nnz_idx] = -A_sparse_val_array[nnz_idx] * G_val.real();
//
//			/* Add identity matrix */
//			if (row_idx==col_idx){
//				IAG_sparse_val_array_cmplx[nnz_idx] += 1;
//			}
//
//			/* Dense test case */
//			IAG_dense_array_cmplx[row_idx*dense_dim + col_idx] = IAG_sparse_val_array_cmplx[nnz_idx];
//		}
//	}
//
//	MKL_INT rows [dense_dim+1];
//	for (int row_idx=0; row_idx<dense_dim+1; row_idx++){
//		rows[row_idx] = A_idx_row_array_csr[row_idx];
//	}
//	MKL_INT cols [A_sparse_dim];
//	for (int col_idx=0; col_idx<A_sparse_dim; col_idx++){
//		cols[col_idx] = A_sparse_col_array[col_idx];
//	}
//
//	MKL_INT sparse_dim_MKL = A_sparse_dim;
//	MKL_INT dense_dim_MKL = dense_dim;
//
//	//double sols_U_col [dense_dim];
//	std::complex<double> sols_U_col [dense_dim];
//
//	/* Solve Faddeev using MKL DSS-routines */
//	auto timestamp_1 = std::chrono::system_clock::now();
//	solve_MM_sparse(IAG_sparse_val_array_cmplx,
//                 	rows,//(MKL_INT*) A_idx_row_array_csr,
//                 	cols,//(MKL_INT*) A_sparse_col_array,
//                 	sparse_dim_MKL,//(MKL_INT) A_sparse_dim,
//                 	A_on_shell_col_array,
//                 	dense_dim_MKL,//(MKL_INT) dense_dim,
//					sols_U_col);
//	auto timestamp_2 = std::chrono::system_clock::now();
//	solve_MM(IAG_dense_array_cmplx, A_dense_array, dense_dim);
//	auto timestamp_3 = std::chrono::system_clock::now();
//	std::chrono::duration<double> time_sparse = timestamp_2 - timestamp_1;
//	std::chrono::duration<double> time_dense  = timestamp_3 - timestamp_2;
//	printf("   - Time used sparse: %.6f\n", time_sparse.count());
//	printf("   - Time used dense:  %.6f\n", time_dense.count());
//
//	for (int row=0; row<dense_dim; row++){
//		//double val_sparse = sols_U_col[row];
//		//double val_dense  = IAG_dense_array_cmplx[row*dense_dim + idx_on_shell];
//		std::complex<double> val_sparse = sols_U_col[row];
//		std::complex<double> val_dense  = IAG_dense_array_cmplx[row*dense_dim + idx_on_shell];
//		//if (val_sparse!=val_dense and std::isnan(val_dense.real())!=true and std::isnan(val_sparse.real())!=true
//		//and std::isnan(val_dense.imag())!=true and std::isnan(val_sparse.imag())!=true){
//			std::cout << "Mismatch. Sparse val: " << val_sparse.real() << " " << val_sparse.imag() << std::endl;
//			std::cout << "Mismatch. Dense val:  " << val_dense.real() << " " << val_dense.imag() << std::endl;
//		//	raise_error("Exit");
//		//}
//	}
//
//	/* Delete temporary arrays */
//	delete [] A_sparse_val_array;
//	delete [] A_sparse_row_array;
//	delete [] A_sparse_col_array;
//	delete [] IAG_sparse_val_array_cmplx;
//}
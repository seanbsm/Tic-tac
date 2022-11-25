
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
					/* NOTE THAT P = P123 + P132 = 2*P123 FOR ANTISYMMETRIC PAIR-STATES */
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
	
	//auto timestamp_start = std::chrono::system_clock::now();
	calculate_PVC_col(PVC_col,
					  idx_alpha_c, idx_p_c, idx_q_c,
					  Nalpha,      Nq_WP,   Np_WP,
					  VC_CM_array,
					  P123_val_array,
					  P123_row_array,
					  P123_col_array,
					  P123_dim);
	//auto timestamp_end = std::chrono::system_clock::now();
	//std::chrono::duration<double>  time1 = timestamp_end - timestamp_start;
	//printf("TIME PVC:  %.6f \n", time1.count()); fflush(stdout);

	/* THOUGHT:
	 * MOVE ALPHA_I OUTWARDS AND GO BACK TO DIRECT APPEND TO COL_ARRAY.
	 * BUT USE NON-ZERO IF-TESTING. */

	/* TOUGHT (CONTRADICTORY TO DIRECT APPEND):
	 * Use dense mat-vec multiplication for sub-blocks  (Can let q be columns of right-vectors? Appealing use of MM-multiplication)
	 * Use dense vec-vec (or vec-mat-vec?) multiplication for A_An */

	/* Generate (C^T x PVC)-column */
	double* CT_subarray     = NULL;
	double* CT_subarray_row = NULL;
	double* PVC_subcol 		= NULL;

	///*  Looping based on matrix-vector product */
	//double* CPVC_subcol 	= NULL;
	//for (size_t idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
	//	/* Beginning of inner-product loops (index "i") */
	//	for (size_t idx_alpha_i=0; idx_alpha_i<Nalpha; idx_alpha_i++){
	//		size_t idx_CT_2N_block = idx_alpha_r*Nalpha + idx_alpha_i;
	//		CT_subarray = CT_RM_array[idx_CT_2N_block];
	//		/* Only do inner-product if CT is not zero due to conservation laws */
	//		if (CT_subarray!=NULL){
	//			for (size_t idx_q_r=0; idx_q_r<Nq_WP; idx_q_r++){
	//				PVC_subcol  = &PVC_col[idx_alpha_i*Nq_WP*Np_WP + idx_q_r*Np_WP];
	//				CPVC_subcol = &col_array[idx_alpha_r*Nq_WP*Np_WP + idx_q_r*Np_WP];
	//				/* Do inner-product over momentum p */
	//				dot_MV(CT_subarray, PVC_subcol, CPVC_subcol, Np_WP, Np_WP);
	//			}
	//		}
	//	}
	//}

	///*  Looping based on first nnz-lookup of PVC elements */
	//for (size_t idx_alpha_i=0; idx_alpha_i<Nalpha; idx_alpha_i++){
	//	for (size_t idx_q_r=0; idx_q_r<Nq_WP; idx_q_r++){
	//		for (size_t idx_p_i=0; idx_p_i<Np_WP; idx_p_i++){
	//			double PVC_element = PVC_col[idx_alpha_i*Nq_WP*Np_WP + idx_q_r*Np_WP + idx_p_i];
	//			if (PVC_element!=0){
	//				for (size_t idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
	//					size_t idx_CT_2N_block = idx_alpha_r*Nalpha + idx_alpha_i;
	//					CT_subarray = CT_RM_array[idx_CT_2N_block];
	//					if (CT_subarray!=NULL){
	//						for (size_t idx_p_r=0; idx_p_r<Np_WP; idx_p_r++){
	//							double prod = PVC_element  * CT_subarray[idx_p_r*Np_WP + idx_p_i];
	//							size_t idx_CPVC = idx_alpha_r*Nq_WP*Np_WP + idx_q_r*Np_WP + idx_p_r;
	//							col_array[idx_CPVC] += prod;
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
	//}
	
	/*  Looping based on first nnz-lookup of CT-matrices */
	/* Loop over rows of col_array */
	//timestamp_start = std::chrono::system_clock::now();
	for (auto idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
		/* Beginning of inner-product loops (index "i") */
		for (auto idx_alpha_i=0; idx_alpha_i<Nalpha; idx_alpha_i++){
			auto idx_CT_2N_block = idx_alpha_r*Nalpha + idx_alpha_i;
			CT_subarray = CT_RM_array[idx_CT_2N_block];
			/* Only do inner-product if CT is not zero due to conservation laws */
			if (CT_subarray!=NULL){
				PVC_subcol = &PVC_col[idx_alpha_i*Nq_WP*Np_WP];
				for (auto idx_q_r=0; idx_q_r<Nq_WP; idx_q_r++){
					for (auto idx_p_i=0; idx_p_i<Np_WP; idx_p_i++){
						//size_t idx_PVC     = idx_alpha_i*Nq_WP*Np_WP + idx_q_r*Np_WP + idx_p_i;
						auto PVC_element = PVC_subcol[idx_q_r*Np_WP + idx_p_i];
						if (PVC_element!=0){
							for (auto idx_p_r=0; idx_p_r<Np_WP; idx_p_r++){
								/* I'm not sure if this is the fastest ordering of the loops */
								//double CT_element  = CT_subarray[idx_p_r*Np_WP + idx_p_i];
								//auto prod = PVC_element  * CT_subarray[idx_p_r*Np_WP + idx_p_i];
								//if (prod!=0){
									auto idx_CPVC = idx_alpha_r*Nq_WP*Np_WP + idx_q_r*Np_WP + idx_p_r;
									col_array[idx_CPVC] += PVC_element  * CT_subarray[idx_p_r*Np_WP + idx_p_i];
									//int nnz_idx = row_to_nnz_array[idx_CPVC];
									//if (nnz_idx==-1){
									//	row_to_nnz_array[idx_CPVC] = num_nnz;
									//	nnz_to_row_array[num_nnz]  = idx_CPVC;
									//	nnz_idx = num_nnz;
									//	num_nnz += 1;
									//}
								//}
							}
						}
					}
				}
			}
		}
	}
	
	//timestamp_end = std::chrono::system_clock::now();
	//std::chrono::duration<double>  time2 = timestamp_end - timestamp_start;
	//printf("TIME CPVC:  %.6f \n", time2.count()); fflush(stdout);
	
	delete [] PVC_col;
}

void calculate_all_CPVC_rows(double*  row_arrays,
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

		/*for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
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
		}*/

		for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
			for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){

				size_t idx_q_com		  = j;

				/* Nucleon-deuteron on-shell (NDOS) indices
				 * (deuteron bound-state p-index is alwasy 0 due to eigenvalue ordering in SWP construction) */
				size_t idx_alpha_NDOS_row = deuteron_idx_array[idx_d_row];
				size_t idx_alpha_NDOS_col = deuteron_idx_array[idx_d_col];
				size_t idx_p_NDOS 	  	  = 0;
				size_t idx_q_NDOS 	   	  = q_com_idx_array[idx_q_com];

				size_t idx_row_NDOS = idx_alpha_NDOS_row*Nq_WP*Np_WP + idx_q_NDOS*Np_WP + idx_p_NDOS;
				size_t idx_col_NDOS = idx_alpha_NDOS_col*Nq_WP*Np_WP + idx_q_NDOS*Np_WP + idx_p_NDOS;

				cdouble U_val = R_array[idx_row_NDOS*dense_dim + idx_col_NDOS];

				size_t idx_NDOS = idx_d_row*num_deuteron_states*num_q_com + idx_d_col*num_q_com + idx_q_com;

				/* Set U-matrix element equal "best" PA */
				U_array[idx_NDOS] = U_val;

				printf("   - U-matrix element for alpha'=%d, alpha=%d, q=%d: %.10e + %.10ei \n", idx_alpha_NDOS_row, idx_alpha_NDOS_col, idx_q_NDOS, U_array[idx_NDOS].real(), U_array[idx_NDOS].imag());
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
					   cdouble*  U_BU_array,
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
					   size_t    P123_sparse_dim,
					   channel_os_indexing chn_os_indexing,
					   run_params run_parameters,
					   std::string file_identification){
						   
	/* Print Pade-approximant convergences */
	bool print_PA_convergences = false;
	/* Print Neumann terms */
	bool print_neumann_terms   = false;
	/* Store Neumann terms */
	bool store_neumann_terms   = true;
	/* Store An-matrices */
	bool store_An_arrays 	   = false;
	/* Store CPVC-kernel to disk */
	bool store_CPVC_array	   = false;
	/* Store CPVC-kernel in memory (faster but much more memory intensive) */
	bool keep_CPVC_in_mem	   = false;
	
	/* Timekeeping variables */
	auto timestamp_start = std::chrono::system_clock::now();
	auto timestamp_end   = std::chrono::system_clock::now();

	/* Use as many threads as possible in MKL-GEMM */
	mkl_set_num_threads(omp_get_max_threads());

	//MKL_NUM_THREADS = omp_get_max_threads();
	//printf("   - Will run with %d MKL-threads \n",mkl_get_max_threads()); fflush(stdout);

	/* Number of on-shell nucleon-deuteron channels (deuteron states can mix, hence ^2) */
	size_t num_on_shell_A_rows = num_deuteron_states * num_q_com;
	size_t num_EL_A_vals 	   = num_deuteron_states * num_deuteron_states * num_q_com;	// Elastic elements
	size_t num_BU_A_vals 	   = num_deuteron_states * chn_os_indexing.num_BU_chns;	// Breakup elements
	
	/* Upper limit on polynomial approximation of Faddeev eq. */
	size_t NM_max = 14;
	size_t num_neumann_terms = 2*NM_max+1;

	/* Coefficients for calculating Pade approximant */
	cdouble* a_coeff_array 	  = new cdouble [ num_neumann_terms * num_EL_A_vals];
	cdouble* a_BU_coeff_array = new cdouble [ num_neumann_terms * num_BU_A_vals];

	/* Dense dimension of 3N-channel */
	size_t dense_dim = Nalpha * Nq_WP * Np_WP;

	/* Allocate row-arrays for A*A^n, where A=(C^T)(P)(VC) 
	 * _array: current iteration of Neumann terms
	 * _array_prev: previous iteration of Neumann terms
	 * _array_comp: compact, past iteration of Neumann terms (compactified to contain only non-converged on-shell elements for faster gemm)
	 * _array_prod: compact, next iteration of Neumann terms (compactified to contain only non-converged on-shell elements for faster gemm) */
	double* re_A_An_row_array 	   = new double [dense_dim * num_on_shell_A_rows];
	double* im_A_An_row_array 	   = new double [dense_dim * num_on_shell_A_rows];
	double* re_A_An_row_array_prev = new double [dense_dim * num_on_shell_A_rows];
	double* im_A_An_row_array_prev = new double [dense_dim * num_on_shell_A_rows];
	double* re_A_An_row_array_comp = new double [dense_dim * num_on_shell_A_rows];
	double* im_A_An_row_array_comp = new double [dense_dim * num_on_shell_A_rows];
	double* re_A_An_row_array_prod = new double [dense_dim * num_on_shell_A_rows];
	double* im_A_An_row_array_prod = new double [dense_dim * num_on_shell_A_rows];
	
	/* Set A_An-arrays to zero */
	for (size_t i=0; i<dense_dim*num_on_shell_A_rows; i++){
		re_A_An_row_array[i] 	  = 0;
		im_A_An_row_array[i] 	  = 0;
		re_A_An_row_array_prev[i] = 0;
		im_A_An_row_array_prev[i] = 0;
		re_A_An_row_array_comp[i] = 0;
		im_A_An_row_array_comp[i] = 0;
		re_A_An_row_array_prod[i] = 0;
		im_A_An_row_array_prod[i] = 0;
	}

	/* Mapping vector from non-converged, compactified rows to full row-storage */
	size_t* A_An_indexing_array = new size_t [num_on_shell_A_rows];

	/* File-paths for storing A_An_row_array and on-shell neumann terms */
	std::string A_An_row_filename      = run_parameters.output_folder + "/An_rows" + file_identification + ".txt";
	std::string neumann_terms_filename = run_parameters.output_folder + "/neumann_terms" + file_identification + ".txt";

	/* Arrays to store Pade-approximants (PA) for each on-shell elastic elements */
	cdouble* pade_approximants_array      = new cdouble [num_EL_A_vals * (NM_max+1)];
	size_t*  pade_approximants_idx_array  = new size_t  [num_EL_A_vals];
	bool*    pade_approximants_conv_array = new bool    [num_EL_A_vals];
	size_t	 num_converged_elements		  = 0;

	/* Arrays to store Pade-approximants (PA) for each on-shell breakup elements */
	cdouble* pade_approximants_BU_array      = new cdouble [num_BU_A_vals * (NM_max+1)];
	size_t*  pade_approximants_BU_idx_array  = new size_t  [num_BU_A_vals];
	bool*    pade_approximants_BU_conv_array = new bool    [num_BU_A_vals];
	size_t	 num_converged_BU_elements		 = 0;

	/* Define CPVC-chunks size */
	size_t num_Gbytes_per_chunk  = 4;
	size_t num_bytes_per_chunk   = num_Gbytes_per_chunk * std::pow(1024,3);
    size_t num_bytes_in_CPVC_col = (sizeof(cdouble) * dense_dim);
    size_t num_cols_per_chunk    = num_bytes_per_chunk / num_bytes_in_CPVC_col;
    size_t num_chunks            = dense_dim / num_cols_per_chunk + 1;
    size_t block_size 			 = num_cols_per_chunk;

	//printf("   - Will calculate kernels in %d chunks of %d GB \n", num_chunks, num_Gbytes_per_chunk); fflush(stdout);

	/* From test-script to program notation */
	size_t    max_num_cols_in_mem = block_size;
	size_t	  num_col_chunks	  = num_chunks;
	double*   CPVC_cols_array     = new double [dense_dim * max_num_cols_in_mem];
	
	/* Allocate row- and column-arrays for (C^T)(P)(VC) */
	int 	 num_threads			    = omp_get_max_threads();
	double*  omp_CPVC_col_array  		= new double [dense_dim * num_threads];
	int*     omp_CPVC_row_to_nnz_array  = new int    [dense_dim * num_threads];
	int*     omp_CPVC_nnz_to_row_array  = new int    [dense_dim * num_threads];

	for (size_t idx_NDOS=0; idx_NDOS<num_EL_A_vals; idx_NDOS++){
		pade_approximants_conv_array[idx_NDOS] = false;
	}
	for (size_t idx_NDOS=0; idx_NDOS<num_BU_A_vals; idx_NDOS++){
		pade_approximants_BU_conv_array[idx_NDOS] = false;
	}

	/* CPVC sparse arrays, used if keep_CPVC_in_mem=true */
	double* CPVC_v_array   = NULL;
	int*    CPVC_c_array   = NULL;
	int*    CPVC_r_array   = NULL;
	size_t* CPVC_csc_array = NULL;
	long long int* CPVC_c_array_LL   = NULL;
	//long long int* CPVC_r_array_LL   = NULL;
	long long int* CPVC_csc_array_LL = NULL;

	
	/* Precalculate kernel CPVC and keep in computer memory, if the option is selected
	 * !!! WARNING: EXTREMELY MEMORY INTENSIVE !!! */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	if (keep_CPVC_in_mem){
		printf("     - Precalculating CPVC-kernel in CSR-sparse format \n"); fflush(stdout);

		size_t CPVC_dim_est = P123_sparse_dim * 100;
		double CPVC_GB_est  = (double) CPVC_dim_est * (2*sizeof(int) + sizeof(double)) / std::pow(1024,3);
		printf("       - Estimated kernel size: %zu non-zero elements (%.2f GB in COO format)\n", CPVC_dim_est, CPVC_GB_est); fflush(stdout);

		/* CPVC-sparse arrays in COO format PER THREAD */
		size_t*  omp_CPVC_dim 	  = new size_t  [num_threads];
		size_t*  omp_CPVC_nnz 	  = new size_t  [num_threads];
		double** omp_CPVC_v_array = new double* [num_threads];
		int**    omp_CPVC_c_array = new int*    [num_threads];
		int**    omp_CPVC_r_array = new int*    [num_threads];
		for (int i=0; i<num_threads; i++){
			size_t omp_dim_est  = P123_sparse_dim * 10 / num_threads; // This is an estimate based on what we see typically, distribute among threads
			omp_CPVC_dim[i]     = omp_dim_est;
			omp_CPVC_nnz[i]     = 0;
			omp_CPVC_v_array[i] = new double [omp_dim_est];
			omp_CPVC_c_array[i] = new int    [omp_dim_est];
			omp_CPVC_r_array[i] = new int    [omp_dim_est];
		}

		/* Fill omp-arrays with nnz elements of CPVC-kernel */
		printf("       - Precalculating ... \n"); fflush(stdout);
		timestamp_start = std::chrono::system_clock::now();
		#pragma omp parallel
		{
			size_t  thread_idx        = omp_get_thread_num();
			double* col_array  	   	  = new double [dense_dim];
			int*    row_to_nnz_array  = new int    [dense_dim];
			int*    nnz_to_row_array  = new int    [dense_dim];

			/* Manually set up parallel for loop, this lets me avoid having to sort columns after matrix construction */
			size_t idx_col_start = (dense_dim/num_threads) *  thread_idx;
			size_t idx_col_end   = (dense_dim/num_threads) * (thread_idx+1) + (thread_idx==num_threads-1)*(dense_dim%num_threads);

			for (size_t idx_col=idx_col_start; idx_col<idx_col_end; idx_col++){
				for (size_t i=0; i<dense_dim; i++){
					col_array[i] = 0;
					row_to_nnz_array[i] = -1;
					nnz_to_row_array[i] = -1;
				}

				size_t idx_alpha_c =  idx_col / (Np_WP*Nq_WP);
				size_t idx_q_c     = (idx_col % (Np_WP*Nq_WP)) /  Np_WP;
				size_t idx_p_c     =  idx_col %  Np_WP;

				/* Calculate CPVC-column */
				size_t num_nnz = 0;
				calculate_CPVC_col(col_array,
								   row_to_nnz_array,
								   nnz_to_row_array,
								   num_nnz,
								   idx_alpha_c, idx_p_c, idx_q_c,
								   Nalpha, Nq_WP, Np_WP,
								   CT_RM_array,
								   VC_CM_array,
								   P123_sparse_val_array,
								   P123_sparse_row_array,
								   P123_sparse_col_array,
								   P123_sparse_dim);


				//std::cout << thread_idx << " " << 1 << std::endl;
				
				/* Lengthen array by 100*dense_dim if the array cannot fit another dense_dim nnz-elements (i.e. max nnz elements from next loop-iteration) */
				if ( omp_CPVC_nnz[thread_idx]+num_nnz+dense_dim >= omp_CPVC_dim[thread_idx] ){
					size_t steplength = 100*dense_dim;
					increase_sparse_array_size(&omp_CPVC_v_array[thread_idx], omp_CPVC_dim[thread_idx], steplength);
					increase_sparse_array_size(&omp_CPVC_r_array[thread_idx], omp_CPVC_dim[thread_idx], steplength);
					increase_sparse_array_size(&omp_CPVC_c_array[thread_idx], omp_CPVC_dim[thread_idx], steplength);
					omp_CPVC_dim[thread_idx] += steplength;
				}

				//std::cout << thread_idx << " " << 2 << std::endl;

				/* Append nnz elements to sparse omp-array */
				size_t curr_nnz = omp_CPVC_nnz[thread_idx];
				for (int nnz=0; nnz<num_nnz; nnz++){
					//std::cout << nnz << " " << curr_nnz << std::endl;
					//std::cout << col_array[nnz_to_row_array[nnz]] << std::endl;
					//std::cout << nnz_to_row_array[nnz] << std::endl;
					//std::cout << idx_col << std::endl;
					omp_CPVC_v_array[thread_idx][curr_nnz+nnz] = col_array[nnz_to_row_array[nnz]];
					omp_CPVC_r_array[thread_idx][curr_nnz+nnz] = nnz_to_row_array[nnz];
					omp_CPVC_c_array[thread_idx][curr_nnz+nnz] = idx_col;
				}
				omp_CPVC_nnz[thread_idx] += num_nnz;
				//std::cout << thread_idx << " " << 3 << std::endl;
			}
		}
		size_t CPVC_num_nnz = 0;
		for (int i=0; i<num_threads; i++){
			CPVC_num_nnz += omp_CPVC_nnz[i];
		}
		double CPVC_GB_true  = (double) CPVC_num_nnz * (2*sizeof(int) + sizeof(double)) / std::pow(1024,3);
		double CPVC_density  = (double) 100.0 * CPVC_num_nnz / std::pow(dense_dim, 2);
		printf("         - Actual kernel size: %zu non-zero elements (%.2f GB in COO format) (%.2f%% density)\n", CPVC_num_nnz, CPVC_GB_true, CPVC_density); fflush(stdout);
		timestamp_end = std::chrono::system_clock::now();
		std::chrono::duration<double>  time_CPVC_construction = timestamp_end - timestamp_start;
		printf("         - Time spent:     %.6f \n", time_CPVC_construction.count()); fflush(stdout);
		printf("         - Done. \n"); fflush(stdout);

		/* Consolidate omp-arrays into a single COO-format array
		 * NOTE: This can be divided into three sections that deallocate e.g. omp-value array
		 * before allocating col-array such that memory is used more efficiently. */
		printf("       - Consolidating distributed arrays ... \n"); fflush(stdout);
		timestamp_start = std::chrono::system_clock::now();

		/* Consolidate values and deallocate parallel memory */
		CPVC_v_array = new double [CPVC_num_nnz];
		size_t idx_nnz = 0;
		for (int i=0; i<num_threads; i++){
			for (size_t j=0; j<omp_CPVC_nnz[i]; j++){
				CPVC_v_array[idx_nnz] = omp_CPVC_v_array[i][j];
				idx_nnz += 1;
			}
		}
		for (int i=0; i<num_threads; i++){
			delete [] omp_CPVC_v_array[i];
		}

		/* Consolidate column-indices and deallocate parallel memory  */
		CPVC_c_array = new int [CPVC_num_nnz];
		idx_nnz = 0;
		for (int i=0; i<num_threads; i++){
			for (size_t j=0; j<omp_CPVC_nnz[i]; j++){
				CPVC_c_array[idx_nnz] = omp_CPVC_c_array[i][j];
				idx_nnz += 1;
			}
		}
		for (int i=0; i<num_threads; i++){
			delete [] omp_CPVC_c_array[i];
		}

		/* Consolidate row-indices and deallocate parallel memory  */
		CPVC_r_array = new int [CPVC_num_nnz];
		idx_nnz = 0;
		for (int i=0; i<num_threads; i++){
			for (size_t j=0; j<omp_CPVC_nnz[i]; j++){
				CPVC_r_array[idx_nnz] = omp_CPVC_r_array[i][j];
				idx_nnz += 1;
			}
		}
		for (int i=0; i<num_threads; i++){
			delete [] omp_CPVC_r_array[i];
		}

		/* Deallocate parallel pointer-arrays */
		delete [] omp_CPVC_dim;
		delete [] omp_CPVC_nnz;
		delete [] omp_CPVC_v_array;
		delete [] omp_CPVC_c_array;
		delete [] omp_CPVC_r_array;

		timestamp_end = std::chrono::system_clock::now();
		std::chrono::duration<double>  time_CPVC_consolidation = timestamp_end - timestamp_start;
		printf("         - Time spent:     %.6f \n", time_CPVC_consolidation.count()); fflush(stdout);
		printf("         - Done. \n"); fflush(stdout);

		///* Sort CPVC COO-sparse array */
		//printf("     - Sorting COO-entries ... \n"); fflush(stdout);
		//timestamp_start = std::chrono::system_clock::now();
		//unsorted_sparse_to_coo_row_major_sorter(&CPVC_v_array,
		//									 	&CPVC_r_array,
		//									 	&CPVC_c_array,
		//									 	CPVC_num_nnz,
		//									 	dense_dim);
		//timestamp_end = std::chrono::system_clock::now();
		//std::chrono::duration<double>  time_CPVC_sorting = timestamp_end - timestamp_start;
		//printf("       - Time spent:     %.6f \n", time_CPVC_sorting.count()); fflush(stdout);
		//printf("       - Done. \n"); fflush(stdout);

		///* Convert COO-array to CSC-array */
		//printf("     - Converting from COO to CSC ... \n"); fflush(stdout);
		//timestamp_start = std::chrono::system_clock::now();
		//CPVC_csc_array = new size_t [dense_dim + 1];
		//coo_to_csc_format_converter(CPVC_r_array,
		//						 	  CPVC_csc_array,
		//						 	  CPVC_num_nnz,
		//						 	  dense_dim);
		//delete [] CPVC_r_array;
		//timestamp_end = std::chrono::system_clock::now();
		//std::chrono::duration<double>  time_CPVC_convertion = timestamp_end - timestamp_start;
		//printf("       - Time spent:     %.6f \n", time_CPVC_convertion.count()); fflush(stdout);
		//printf("       - Done. \n"); fflush(stdout);
		///* Copy column-array from int to long long int */
		//CPVC_c_array_LL   = new long long int [CPVC_num_nnz];
		//for (size_t i=0; i<CPVC_num_nnz; i++){
		//	CPVC_c_array_LL[i] = CPVC_c_array[i];
		//}
		//delete [] CPVC_c_array;
		///* Copy CSC-array from int to long long int */
		//CPVC_csc_array_LL = new long long int [dense_dim+1];
		//for (int i=0; i<dense_dim+1; i++){
		//	CPVC_csc_array_LL[i] = CPVC_csc_array[i];
		//}
		//delete [] CPVC_csc_array;

		std::string filename = "kernel.h5";   //P123_folder + "/CPVC_"
						 				   //+ to_string(two_J_3N) + "_" + to_string(P_3N)
						 				   //+ "_Np_" + to_string(Np_WP) + "_Nq_" + to_string(Nq_WP)
						 				   //+ "_J2max_" + to_string(J_2N_max) + "_TFC_" + to_string(thread_idx)
						 				   //+ "_" + to_string(current_TFC) + ".h5";

		/* Store sparse CPVC-array to file */
		printf("       - Writing kernel to h5 ... \n"); fflush(stdout);
		timestamp_start = std::chrono::system_clock::now();
		store_sparse_matrix_h5(CPVC_v_array,
							   CPVC_r_array,
							   CPVC_c_array,
							   CPVC_num_nnz,
							   dense_dim,
					   		   filename,
							   false);
		timestamp_end = std::chrono::system_clock::now();
		std::chrono::duration<double>  time_CPVC_storage = timestamp_end - timestamp_start;
		printf("         - Successs. \n"); fflush(stdout);
		printf("         - Time spent:     %.6f \n", time_CPVC_storage.count()); fflush(stdout);
		printf("         - Done. \n"); fflush(stdout);

		/* Read sparse CPVC-array from file */
		printf("       - Reading kernel from h5 ... \n"); fflush(stdout);
		timestamp_start = std::chrono::system_clock::now();
		delete [] CPVC_v_array;
		CPVC_v_array = NULL;
		delete [] CPVC_r_array;
		CPVC_r_array = NULL;
		delete [] CPVC_c_array;
		CPVC_c_array = NULL;
		read_sparse_matrix_h5(&CPVC_v_array,
							   &CPVC_r_array,
							   &CPVC_c_array,
							   CPVC_num_nnz,
							   dense_dim,
					   		   filename,
							   false);
		timestamp_end = std::chrono::system_clock::now();
		std::chrono::duration<double>  time_CPVC_reading = timestamp_end - timestamp_start;
		printf("         - Successs. \n"); fflush(stdout);
		printf("         - Time spent:     %.6f \n", time_CPVC_reading.count()); fflush(stdout);
		printf("         - Done. \n"); fflush(stdout);

		printf("       - Done. \n"); fflush(stdout);
	}
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* ################################################################################################################### */

	/* Set initial values for A_Kn_row_array, where K^n=1 for n=0 */
	printf("     - Working on Pade approximant P[N,M] for N=%d, M=%d \n",0,0); fflush(stdout);
	printf("       - Calculating on-shell rows of A*K^n for n=%d. \n", 0); fflush(stdout);
	timestamp_start = std::chrono::system_clock::now();
	/* Calculate CPVC-row and write to A_Kn_row_array_prev */
	calculate_all_CPVC_rows(re_A_An_row_array_prev,
							q_com_idx_array, num_q_com,
			   				deuteron_idx_array, num_deuteron_states,
							Nalpha, Nq_WP, Np_WP,
							CT_RM_array,
							VC_CM_array,
							P123_sparse_val_array,
							P123_sparse_row_array,
							P123_sparse_col_array,
							P123_sparse_dim);
	timestamp_end = std::chrono::system_clock::now();
	std::chrono::duration<double> time = timestamp_end - timestamp_start;
	printf("         - Time generating CPVC-rows:     %.6f \n", time.count()); fflush(stdout);
	printf("         - Done \n"); fflush(stdout);
	
	/* First Neumann-term */
	printf("       - Extracting on-shell Neumann-series terms a_n=A*K^n for n=%d. \n",0); fflush(stdout);
	/* Extract elastic terms */
	for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
		for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
			for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
				/* Nucleon-deuteron on-shell (NDOS) indices
				 * (deuteron bound-state p-index is alwasy 0 due to eigenvalue ordering in SWP construction) */
				size_t idx_alpha_NDOS_row = deuteron_idx_array[idx_d_row];
				size_t idx_alpha_NDOS_col = deuteron_idx_array[idx_d_col];
				size_t idx_q_NDOS 	  = q_com_idx_array[idx_q_com];
				size_t idx_p_NDOS 	  = 0;

				size_t idx_row_NDOS   = idx_d_row*num_q_com + idx_q_com;
				size_t idx_col_NDOS   = idx_alpha_NDOS_col*Nq_WP*Np_WP + idx_q_NDOS*Np_WP + idx_p_NDOS;

				/* Calculate coefficient */
				cdouble a_coeff = re_A_An_row_array_prev[idx_row_NDOS*dense_dim + idx_col_NDOS];

				/* Store coefficient */
				size_t idx_NDOS = idx_d_row*num_deuteron_states*num_q_com + idx_d_col*num_q_com + idx_q_com;
				a_coeff_array[idx_NDOS*num_neumann_terms] = a_coeff;
				
				if (print_neumann_terms){
					printf("         - Neumann term %d for alpha'=%d, alpha=%d, q=%d: %.16e + %.16ei \n", 0, idx_alpha_NDOS_row, idx_alpha_NDOS_col, idx_q_NDOS, a_coeff.real(), a_coeff.imag());
					fflush(stdout);
				}
			}
		}
	}
	/* Extract breakup terms */
	for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
		int BU_chn_start = chn_os_indexing.q_com_BU_idx_array[idx_q_com];
		int BU_chn_end   = chn_os_indexing.q_com_BU_idx_array[idx_q_com+1];
		for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
			for (size_t idx_BU_chn=BU_chn_start; idx_BU_chn<BU_chn_end; idx_BU_chn++){
				
				size_t idx_alpha_NDOS = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 0];
				size_t idx_q_NDOS 	  = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 1];
				size_t idx_p_NDOS 	  = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 2];

				size_t idx_row_NDOS   = idx_d_row*num_q_com + idx_q_com;
				size_t idx_col_NDOS   = idx_alpha_NDOS*Nq_WP*Np_WP + idx_q_NDOS*Np_WP + idx_p_NDOS;

				/* Calculate coefficient */
				cdouble a_BU_coeff = re_A_An_row_array_prev[idx_row_NDOS*dense_dim + idx_col_NDOS];
				
				/* Store coefficient */
				size_t idx_NDOS = idx_d_row*chn_os_indexing.num_BU_chns + idx_BU_chn;
				a_BU_coeff_array[idx_NDOS*num_neumann_terms] = a_BU_coeff;
			}
		}
	}
	printf("         - Done \n"); fflush(stdout);

	if (store_An_arrays){
		printf("       - Storing matrix A*K^n for n=%d to output-folder. \n", 0); fflush(stdout);
		std::string array_seperator_text = "n = " + std::to_string(0);
		store_sep_complex_matrix(re_A_An_row_array_prev,
								 im_A_An_row_array_prev,
    	                         num_on_shell_A_rows,
						         dense_dim,
						         dense_dim,
						         A_An_row_filename,
						         true,
							     array_seperator_text);
		printf("         - Done \n");
	}
	if (store_neumann_terms){
		printf("       - Storing on-shell Neumann-series terms a_n=A*K^n for n=%d to output-folder. \n", 0); fflush(stdout);
		std::string array_seperator_text = "n = " + std::to_string(0);
		store_complex_matrix(a_coeff_array,
    	                     num_deuteron_states*num_deuteron_states*num_q_com,
			    			 num_neumann_terms,
			    			 num_neumann_terms,
			    			 neumann_terms_filename,
			    			 true,
							 array_seperator_text);
		printf("         - Done \n");
	}

	/* Loop over number of Pade-terms we use */
	for (size_t NM=0; NM<NM_max+1; NM++){
		if (num_converged_elements==num_EL_A_vals){
			printf("     - Convergence reached for all on-shell elements! \n"); fflush(stdout);
			break;
		}

		if (NM!=0){
			printf("     - Working on Pade approximant P[N,M] for N=%d, M=%d \n",NM,NM); fflush(stdout);
		}
		
		size_t counter_array [100];
		for (size_t i=0; i<100; i++){
			counter_array[i] = 0;
		}

		/* Time-keeper array for parallel environment */
		double*  times_array = new double [3*num_threads];
		
		for (int n=2*NM-1; n<2*NM+1; n++){
			printf("       - Working on Neumann-terms for n=%d. \n", n); fflush(stdout);
			/* We've already done n=0 above */
			if (n<=0){
				continue;
			}
			
			/* Initialise time-profile variables */
			double time_resolvent        = 0;
			double time_CPVC_cols        = 0;
			double time_An_CPVC_multiply = 0;
			double time_neumann          = 0;

			double timestamp_neumann_start = omp_get_wtime();

			/* Calculate all a-coefficients for calculated CPVC-column */
			double timestamp_resolvent_start = omp_get_wtime();
			printf("       - Multiplying in resolvent with An. \n"); fflush(stdout);
			for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
				for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
					size_t idx_row_NDOS = idx_d_row*num_q_com + idx_q_com;

					/* Multiply An by G */
					for (size_t idx=0; idx<dense_dim; idx++){
						double re_An = re_A_An_row_array_prev[idx_row_NDOS*dense_dim + idx];
						double im_An = im_A_An_row_array_prev[idx_row_NDOS*dense_dim + idx];
						double re_G  = G_array[idx_q_com*dense_dim + idx].real();
						double im_G  = G_array[idx_q_com*dense_dim + idx].imag();
						re_A_An_row_array_prev[idx_row_NDOS*dense_dim + idx] = re_An*re_G - im_An*im_G;
						im_A_An_row_array_prev[idx_row_NDOS*dense_dim + idx] = re_An*im_G + im_An*re_G;
					}
				}
			}
			double timestamp_resolvent_end    = omp_get_wtime();
			time_resolvent = timestamp_resolvent_end - timestamp_resolvent_start;

			size_t num_non_conv_rows = 0;
			for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
				for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
					/* Nucleon-deuteron on-shell (NDOS) indices
					 * (deuteron bound-state p-index is alwasy 0 due to eigenvalue ordering in SWP construction) */
					size_t idx_alpha_NDOS_row = deuteron_idx_array[idx_d_row];
					size_t idx_p_NDOS 	  = 0;
					size_t idx_q_NDOS 	  = q_com_idx_array[idx_q_com];

					size_t idx_row_NDOS   = idx_d_row*num_q_com + idx_q_com;

					/* See if row contains unconverged, on-shell elements */
					bool row_conv = true;
					/* Check elastic channels */
					for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
						size_t idx_NDOS = idx_d_row*num_deuteron_states*num_q_com + idx_d_col*num_q_com + idx_q_com;
						if (pade_approximants_conv_array[idx_NDOS]==false){
							row_conv = false;
						}
					}
					/* Check breakup channels */
					int BU_chn_start = chn_os_indexing.q_com_BU_idx_array[idx_q_com];
					int BU_chn_end   = chn_os_indexing.q_com_BU_idx_array[idx_q_com+1];
					for (size_t idx_BU_chn=BU_chn_start; idx_BU_chn<BU_chn_end; idx_BU_chn++){
						size_t idx_NDOS = idx_d_row*chn_os_indexing.num_BU_chns + idx_BU_chn;
						if (pade_approximants_BU_conv_array[idx_NDOS]==false){
							row_conv = false;
						}
					}
					
					if (row_conv==false){
						for (size_t i=0; i<dense_dim; i++){
							re_A_An_row_array_comp[num_non_conv_rows*dense_dim + i] = re_A_An_row_array_prev[idx_row_NDOS*dense_dim + i];
							im_A_An_row_array_comp[num_non_conv_rows*dense_dim + i] = im_A_An_row_array_prev[idx_row_NDOS*dense_dim + i];
						}
						A_An_indexing_array[num_non_conv_rows] = idx_row_NDOS;
						num_non_conv_rows += 1;
					}
				}
			}
			
			printf("       - Calculating on-shell rows of A*K^n for n=%d. \n", n); fflush(stdout);
			if (keep_CPVC_in_mem==false){
				for (size_t idx_col_chunk=0; idx_col_chunk<num_col_chunks; idx_col_chunk++){

					double timestamp_CPVC_chunk_start = omp_get_wtime();

					size_t idx_col_start =  idx_col_chunk    * max_num_cols_in_mem;
					size_t idx_col_end   = (idx_col_chunk+1) * max_num_cols_in_mem;

					if (idx_col_end>dense_dim){
						idx_col_end = dense_dim;
					}

					size_t cols_in_chunk = idx_col_end - idx_col_start;

					/* Reset CPVC-columns when array is filled */
					for (size_t idx=0; idx<dense_dim * max_num_cols_in_mem; idx++){
						CPVC_cols_array[idx] = 0;
					}

					#pragma omp parallel //num_threads(1)
					{

					size_t  thread_idx             = omp_get_thread_num();
					double* CPVC_col_array  	   = &omp_CPVC_col_array	    [thread_idx*dense_dim];
					int*    CPVC_row_to_nnz_array  = &omp_CPVC_row_to_nnz_array [thread_idx*dense_dim];
					int*    CPVC_nnz_to_row_array  = &omp_CPVC_nnz_to_row_array [thread_idx*dense_dim];

					#pragma omp for
					for (size_t idx_col=idx_col_start; idx_col<idx_col_end; idx_col++){
						size_t idx_alpha_c = idx_col / (Np_WP*Nq_WP);
						size_t idx_q_c     = (idx_col % (Np_WP*Nq_WP)) /  Np_WP;
						size_t idx_p_c     = idx_col %  Np_WP;

						/* Calculate CPVC-column */
						size_t CPVC_num_nnz = 0;
						size_t CPVC_col_idx = idx_col - idx_col_start;
						calculate_CPVC_col(&CPVC_cols_array[CPVC_col_idx*dense_dim],//CPVC_col_array,
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
					}
					}
					double timestamp_CPVC_chunk_end   = omp_get_wtime();
					time_CPVC_cols += timestamp_CPVC_chunk_end - timestamp_CPVC_chunk_start;

					double beta  = 0;
					double alpha = 1;
					MKL_INT M    = num_non_conv_rows;// num_on_shell_A_rows
					MKL_INT N    = cols_in_chunk;// max_num_cols_in_mem;
					MKL_INT K    = dense_dim;
					MKL_INT lda  = dense_dim;
					MKL_INT ldb  = dense_dim;//max_num_cols_in_mem;
					MKL_INT ldc  = dense_dim;
					double* re_A = &re_A_An_row_array_comp[0];
					double* im_A = &im_A_An_row_array_comp[0];
					double* B 	 = &CPVC_cols_array[0];
					double* re_C = &re_A_An_row_array_prod[idx_col_start];
					double* im_C = &im_A_An_row_array_prod[idx_col_start];

					double timestamp_gemm_start = omp_get_wtime();
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, alpha, re_A, lda, B, ldb, beta, re_C, ldc);	// real multiplication
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, alpha, im_A, lda, B, ldb, beta, im_C, ldc);	// imag multiplication
					double timestamp_gemm_end   = omp_get_wtime();
					time_An_CPVC_multiply += timestamp_gemm_end - timestamp_gemm_start;
				}
			}
			else{
				double timestamp_gemm_start = omp_get_wtime();
				const char   ordering = 'R';
				const char   trans 	  = 'T';
				const double alpha 	  = 1.0;
				/* Transpose An before sparse multiplication */
				mkl_dimatcopy(ordering, trans, num_non_conv_rows, dense_dim, alpha, re_A_An_row_array_comp, dense_dim, num_non_conv_rows);
				mkl_dimatcopy(ordering, trans, num_non_conv_rows, dense_dim, alpha, im_A_An_row_array_comp, dense_dim, num_non_conv_rows);
				/* Multiply CPVC with An using sparse multiplication */
				//dot_MM_sparse(CPVC_v_array, CPVC_c_array_LL, CPVC_csc_array_LL, re_A_An_row_array_comp, re_A_An_row_array_prod, dense_dim, dense_dim, num_non_conv_rows, true);
				//dot_MM_sparse(CPVC_v_array, CPVC_c_array_LL, CPVC_csc_array_LL, re_A_An_row_array_comp, im_A_An_row_array_prod, dense_dim, dense_dim, num_non_conv_rows, true);
				/* Transpose An+1 after sparse multiplication */
				mkl_dimatcopy(ordering, trans, dense_dim, num_non_conv_rows, alpha, re_A_An_row_array_prod, num_non_conv_rows, dense_dim);
				mkl_dimatcopy(ordering, trans, dense_dim, num_non_conv_rows, alpha, im_A_An_row_array_prod, num_non_conv_rows, dense_dim);
				double timestamp_gemm_end   = omp_get_wtime();
				time_An_CPVC_multiply += timestamp_gemm_end - timestamp_gemm_start;
			}

			/* Write compact format back to full format */
			for (size_t r=0; r<num_non_conv_rows; r++){
				size_t idx_row_NDOS = A_An_indexing_array[r];
				for (size_t i=0; i<dense_dim; i++){
					re_A_An_row_array[idx_row_NDOS*dense_dim + i] = re_A_An_row_array_prod[r*dense_dim + i];
					im_A_An_row_array[idx_row_NDOS*dense_dim + i] = im_A_An_row_array_prod[r*dense_dim + i];
				}
			}
			double timestamp_neumann_end = omp_get_wtime();
			time_neumann = timestamp_neumann_end - timestamp_neumann_start;

			printf("         - Time multiplying An with G:    %.6f s \n", time_resolvent);
			printf("         - Time generating CPVC-cols:     %.6f s \n", time_CPVC_cols);
			printf("         - Time multiplying An with CPVC: %.6f s \n", time_An_CPVC_multiply);
			printf("         - Total time:                    %.6f s \n", time_neumann);
			printf("         - Done \n"); fflush(stdout);

			//size_t nnz_counts = 0;
			//for (size_t i=0; i<100; i++){
			//	nnz_counts += counter_array[i];
			//}
			//printf("NUMBER OF NNZ ELEMENTS IN A: %zu \n", nnz_counts);
			//printf("NUMBER OF NNZ ELEMENTS IN P: %zu \n", P123_sparse_dim);

			/* Rewrite previous A_An with current A_An */
			for (size_t i=0; i<num_on_shell_A_rows*dense_dim; i++){
				re_A_An_row_array_prev[i] = re_A_An_row_array[i];
				im_A_An_row_array_prev[i] = im_A_An_row_array[i];
			}

			/* Extract coefficients "a" for Pade approximant */
			printf("       - Extracting on-shell Neumann-series terms a_n=A*K^n for n=%d. \n", n); fflush(stdout);
			/* Extract elastic terms */
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

						size_t idx_NDOS = idx_d_row*num_deuteron_states*num_q_com + idx_d_col*num_q_com + idx_q_com;

						/* Check if we've already reached convergence for this on-shell element */
						if (pade_approximants_conv_array[idx_NDOS]==true){
							continue;
						}

						/* Calculate coefficient */
						cdouble a_coeff = {re_A_An_row_array[idx_row_NDOS*dense_dim + idx_col_NDOS], im_A_An_row_array[idx_row_NDOS*dense_dim + idx_col_NDOS]};

						/* Store coefficient */
						a_coeff_array[idx_NDOS*num_neumann_terms + n] = a_coeff;

						if (print_neumann_terms){
							printf("         - Neumann term %d for alpha'=%d, alpha=%d, q=%d: %.16e + %.16ei \n", n, idx_alpha_NDOS_row, idx_alpha_NDOS_col, idx_q_NDOS, a_coeff.real(), a_coeff.imag());
							//printf("%d\n", idx_row_NDOS*dense_dim + idx_col_NDOS);
						}
					}
				}
			}
			/* Extract breakup terms */
			for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
				int BU_chn_start = chn_os_indexing.q_com_BU_idx_array[idx_q_com];
				int BU_chn_end   = chn_os_indexing.q_com_BU_idx_array[idx_q_com+1];
				for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
					for (size_t idx_BU_chn=BU_chn_start; idx_BU_chn<BU_chn_end; idx_BU_chn++){
						
						size_t idx_alpha_NDOS = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 0];
						size_t idx_q_NDOS 	  = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 1];
						size_t idx_p_NDOS 	  = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 2];

						size_t idx_row_NDOS   = idx_d_row*num_q_com + idx_q_com;
						size_t idx_col_NDOS   = idx_alpha_NDOS*Nq_WP*Np_WP + idx_q_NDOS*Np_WP + idx_p_NDOS;

						size_t idx_NDOS = idx_d_row*chn_os_indexing.num_BU_chns + idx_BU_chn;

						/* Check if we've already reached convergence for this on-shell element */
						if (pade_approximants_BU_conv_array[idx_NDOS]==true){
							continue;
						}

						/* Calculate coefficient */
						cdouble a_BU_coeff = {re_A_An_row_array[idx_row_NDOS*dense_dim + idx_col_NDOS], im_A_An_row_array[idx_row_NDOS*dense_dim + idx_col_NDOS]};
						
						/* Store coefficient */
						a_BU_coeff_array[idx_NDOS*num_neumann_terms + n] = a_BU_coeff;
					}
				}
			}
			printf("         - Done \n");

			if (store_An_arrays){
			printf("       - Storing matrix A*K^n for n=%d to output-folder. \n", n); fflush(stdout);
			std::string array_seperator_text = "n = " + std::to_string(n);
		    store_sep_complex_matrix(re_A_An_row_array,
									 im_A_An_row_array,
                                     num_on_shell_A_rows,
		    				         dense_dim,
		    				         dense_dim,
		    				         A_An_row_filename,
		    				         false,
								     array_seperator_text);
			printf("         - Done \n");
			}
			if (store_neumann_terms){
				printf("       - Storing on-shell Neumann-series terms a_n=A*K^n for n=%d to output-folder. \n", n); fflush(stdout);
				std::string array_seperator_text = "n = " + std::to_string(n);
				store_complex_matrix(a_coeff_array,
            	                     num_deuteron_states*num_deuteron_states*num_q_com,
		    					     num_neumann_terms,
		    					     num_neumann_terms,
		    					     neumann_terms_filename,
		    					     false,
									 array_seperator_text);
				printf("         - Done \n");
			}
		}
		delete [] times_array;

		printf("       - Calculating Pade approximants PA[%d,%d]. \n", NM, NM); fflush(stdout);
		/* Calculate Pade approximants (PA) for elastic amplitudes */
		for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
			for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
				for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){

					size_t idx_NDOS = idx_d_row*num_deuteron_states*num_q_com + idx_d_col*num_q_com + idx_q_com;

					/* Check and skip if we've already reached convergence for this on-shell element */
					if (pade_approximants_conv_array[idx_NDOS]==true){
						continue;
					}

					/* Calculate and append PA */
					cdouble PA = pade_approximant(&a_coeff_array[idx_NDOS*num_neumann_terms], NM, NM, 1);

					pade_approximants_array[idx_NDOS*(NM_max+1) + NM] = PA;
					
					/* See if we've reached convergence with this iteration */
					size_t idx_best_PA = 0;
					bool convergence_reached = false;

					/* Find minimum PA from previous calculations */
					double min_PA_diff = 1;
					for (int NM_prev=0; NM_prev<NM; NM_prev++){
						cdouble PA_prev = pade_approximants_array[idx_NDOS*(NM_max+1) + NM_prev];

						/* Calculate difference between PAs from previous PA-calculations */
						double PA_diff_prev = std::abs(PA_prev - pade_approximants_array[idx_NDOS*(NM_max+1) + NM_prev-1]);

						/* Ignore PA_diff_prev if numerically equal to the previous PA_diff, overwrite if smaller than min_PA_diff */
						if (PA_diff_prev<min_PA_diff && PA_diff_prev>1e-15){
							idx_best_PA = NM_prev;
							min_PA_diff = PA_diff_prev;
						}
					}
					/* See if current PA is better/worse than previous minimum */
					double PA_diff_curr = std::abs(PA - pade_approximants_array[idx_NDOS*(NM_max+1) + NM - 1]);
					
					/* Ignore PA_diff_curr if numerically equal to the previous PA_diff_curr, overwrite if smaller than min_PA_diff */
					if (PA_diff_curr<min_PA_diff && PA_diff_curr>1e-15){
						idx_best_PA = NM;
						min_PA_diff = PA_diff_curr;
					}

					/* Criterias for convergence¨
					 * If any are fulfilled, we set convergence to true for current on-shell element */
					bool convergence_criteria_0 = (NM==NM_max);																				// Cannot go past max NM
					bool convergence_criteria_1 = (NM-idx_best_PA>4);																		// If nothing better is found in the last 3 PAs, we assume we found the best
					bool convergence_criteria_2 = (min_PA_diff<1e-6*std::abs(pade_approximants_array[idx_NDOS*(NM_max+1) + idx_best_PA]));	// If the difference is less than the 4th significant digit, we assume "good enough"
					bool convergence_criteria_3 = (min_PA_diff<1e-7);																		// If we are below single precision resolution, assume convergence
					
					if (convergence_criteria_0 ||
						convergence_criteria_1 ||
						convergence_criteria_2 ||
						convergence_criteria_3){
						pade_approximants_conv_array[idx_NDOS] = true;
						pade_approximants_idx_array[idx_NDOS]  = idx_best_PA;
						num_converged_elements += 1;
					}

					//if (print_PA_convergences){
					//	printf("PA[%d,%d] = %.16e + %.16ei, PA_diff = %.16e \n", NM,NM,PA.real(), PA.imag(), PA_diff);
					//}
				}
			}
		}
		/* Calculate Pade approximants (PA) for breakup amplitudes */
		for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
			int BU_chn_start = chn_os_indexing.q_com_BU_idx_array[idx_q_com];
			int BU_chn_end   = chn_os_indexing.q_com_BU_idx_array[idx_q_com+1];
			for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
				for (size_t idx_BU_chn=BU_chn_start; idx_BU_chn<BU_chn_end; idx_BU_chn++){

					size_t idx_NDOS = idx_d_row*chn_os_indexing.num_BU_chns + idx_BU_chn;

					/* Check and skip if we've already reached convergence for this on-shell element */
					if (pade_approximants_BU_conv_array[idx_NDOS]==true){
						continue;
					}

					/* Calculate and append PA */
					cdouble PA = pade_approximant(&a_BU_coeff_array[idx_NDOS*num_neumann_terms], NM, NM, 1);

					pade_approximants_BU_array[idx_NDOS*(NM_max+1) + NM] = PA;
					
					/* See if we've reached convergence with this iteration */
					size_t idx_best_PA = 0;
					bool convergence_reached = false;

					/* Find minimum PA from previous calculations */
					double min_PA_diff = 1;
					for (int NM_prev=0; NM_prev<NM; NM_prev++){
						cdouble PA_prev = pade_approximants_BU_array[idx_NDOS*(NM_max+1) + NM_prev];

						/* Calculate difference between PAs from previous PA-calculations */
						double PA_diff_prev = std::abs(PA_prev - pade_approximants_BU_array[idx_NDOS*(NM_max+1) + NM_prev-1]);

						/* Ignore PA_diff_prev if numerically equal to the previous PA_diff, overwrite if smaller than min_PA_diff */
						if (PA_diff_prev<min_PA_diff && PA_diff_prev>1e-15){
							idx_best_PA = NM_prev;
							min_PA_diff = PA_diff_prev;
						}
					}
					/* See if current PA is better/worse than previous minimum */
					double PA_diff_curr = std::abs(PA - pade_approximants_BU_array[idx_NDOS*(NM_max+1) + NM - 1]);
					
					/* Ignore PA_diff_curr if numerically equal to the previous PA_diff_curr, overwrite if smaller than min_PA_diff */
					if (PA_diff_curr<min_PA_diff && PA_diff_curr>1e-15){
						idx_best_PA = NM;
						min_PA_diff = PA_diff_curr;
					}

					/* Criterias for convergence¨
					 * If any are fulfilled, we set convergence to true for current on-shell element */
					bool convergence_criteria_0 = (NM==NM_max);																					// Cannot go past max NM
					bool convergence_criteria_1 = (NM-idx_best_PA>4);																			// If nothing better is found in the last 3 PAs, we assume we found the best
					bool convergence_criteria_2 = (min_PA_diff<1e-6*std::abs(pade_approximants_BU_array[idx_NDOS*(NM_max+1) + idx_best_PA]));	// If the difference is less than the 4th significant digit, we assume "good enough"
					bool convergence_criteria_3 = (min_PA_diff<1e-7);																			// If we are below single precision resolution, assume convergence
					
					if (convergence_criteria_0 ||
						convergence_criteria_1 ||
						convergence_criteria_2 ||
						convergence_criteria_3){
						pade_approximants_BU_conv_array[idx_NDOS] = true;
						pade_approximants_BU_idx_array[idx_NDOS]  = idx_best_PA;
						num_converged_elements += 1;
					}

					//if (print_PA_convergences){
					//	printf("PA[%d,%d] = %.16e + %.16ei, PA_diff = %.16e \n", NM,NM,PA.real(), PA.imag(), PA_diff);
					//}
				}
			}
		}
		printf("         - Done \n"); fflush(stdout);
	}

	printf("     - Extracting on-shell U-matrix elements \n"); fflush(stdout);
	/* Set on-shell elastic U-matrix elements equal "best" PA */
	for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
		for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
			for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
				/* Nucleon-deuteron on-shell (NDOS) indices
				 * (deuteron bound-state p-index is alwasy 0 due to eigenvalue ordering in SWP construction) */
				size_t idx_alpha_NDOS_row = deuteron_idx_array[idx_d_row];
				size_t idx_alpha_NDOS_col = deuteron_idx_array[idx_d_col];
				size_t idx_q_NDOS 	   	  = q_com_idx_array[idx_q_com];

				size_t idx_NDOS = idx_d_row*num_deuteron_states*num_q_com + idx_d_col*num_q_com + idx_q_com;

				size_t idx_best_PA = pade_approximants_idx_array[idx_NDOS];

				U_array[idx_NDOS] = pade_approximants_array[idx_NDOS*(NM_max+1) + idx_best_PA];
				printf("       - U-matrix element for alpha'=%d, alpha=%d, q=%d: %.10e + %.10ei \n", idx_alpha_NDOS_row, idx_alpha_NDOS_col, idx_q_NDOS, U_array[idx_NDOS].real(), U_array[idx_NDOS].imag());
			}
		}
	}
	/* Set on-shell breakup U-matrix elements equal "best" PA */
	for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
		int BU_chn_start = chn_os_indexing.q_com_BU_idx_array[idx_q_com];
		int BU_chn_end   = chn_os_indexing.q_com_BU_idx_array[idx_q_com+1];
		for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
			for (size_t idx_BU_chn=BU_chn_start; idx_BU_chn<BU_chn_end; idx_BU_chn++){
				/* Nucleon-deuteron breakup on-shell (NDOS) indices */
				size_t idx_alpha_NDOS = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 0];
				size_t idx_q_NDOS 	  = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 1];
				size_t idx_p_NDOS 	  = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 2];

				size_t idx_NDOS = idx_d_row*chn_os_indexing.num_BU_chns + idx_BU_chn;

				size_t idx_best_PA = pade_approximants_BU_idx_array[idx_NDOS];

				U_BU_array[idx_NDOS] = pade_approximants_BU_array[idx_NDOS*(NM_max+1) + idx_best_PA];
				//printf("       - U-matrix element for alpha'=%d, alpha=%d, q=%d: %.10e + %.10ei \n", idx_alpha_NDOS_row, idx_alpha_NDOS_col, idx_q_NDOS, U_array[idx_NDOS].real(), U_array[idx_NDOS].imag());
			}
		}
	}
	printf("       - Done \n"); fflush(stdout);

	printf("       - Storing on-shell Neumann-series terms a_n=A*K^n for all n to output-folder. \n"); fflush(stdout);
	store_complex_matrix(a_coeff_array,
                         num_deuteron_states*num_deuteron_states*num_q_com,
					     num_neumann_terms,
					     num_neumann_terms,
					     neumann_terms_filename,
					     true,
						 "Neumann terms");
	printf("         - Done \n");

	/* Free allocated working space */
	delete [] a_coeff_array;
	delete [] a_BU_coeff_array;
	delete [] CPVC_cols_array;
	delete [] omp_CPVC_col_array;
	delete [] omp_CPVC_row_to_nnz_array;
	delete [] omp_CPVC_nnz_to_row_array;
	delete [] re_A_An_row_array;
	delete [] im_A_An_row_array;
	delete [] re_A_An_row_array_prev;
	delete [] im_A_An_row_array_prev;
	delete [] re_A_An_row_array_comp;
	delete [] im_A_An_row_array_comp;
	delete [] re_A_An_row_array_prod;
	delete [] im_A_An_row_array_prod;
	delete [] pade_approximants_array;
	delete [] pade_approximants_idx_array;
	delete [] pade_approximants_conv_array;
	delete [] pade_approximants_BU_array;
	delete [] pade_approximants_BU_idx_array;
	delete [] pade_approximants_BU_conv_array;
}

void solve_faddeev_equations(cdouble*  U_array,
							 cdouble*  U_BU_array,
							 cdouble*  G_array,
							 double*   P123_sparse_val_array,
							 int*      P123_sparse_row_array,
							 size_t*   P123_sparse_col_array_csc,
							 size_t    P123_sparse_dim,
							 double*   V_WP_unco_array,
							 double*   V_WP_coup_array,
							 swp_statespace swp_states,
							 channel_os_indexing chn_os_indexing,
							 pw_3N_statespace pw_states,
							 std::string file_identification,
					         run_params run_parameters){
	
	/* Make local pointers & variables for on-shell channel-indexing */
	int*   q_com_idx_array		= chn_os_indexing.q_com_idx_array;
	int*   deuteron_idx_array	= chn_os_indexing.deuteron_idx_array;
	size_t num_q_com			= (size_t) chn_os_indexing.num_T_lab;
	size_t num_deuteron_states	= (size_t) chn_os_indexing.num_deuteron_states;
	/* Make local pointers & variables for SWP-statespace */
	size_t    Nq_WP					= (size_t) swp_states.Nq_WP;
	size_t    Np_WP					= (size_t) swp_states.Np_WP;
	double*   C_WP_unco_array		= swp_states.C_SWP_unco_array;
	double*   C_WP_coup_array		= swp_states.C_SWP_coup_array;
	int  	  num_2N_unco_states	= swp_states.num_2N_unco_states;
	int  	  num_2N_coup_states	= swp_states.num_2N_coup_states;
	/* Make local variable for pw-statespace */
	size_t Nalpha = pw_states.Nalpha;

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
									   num_2N_unco_states,
									   num_2N_coup_states,
									   Np_WP,
									   pw_states,
									   run_parameters);

	/* Create VC-product pointer-arrays in column-major format */
	double** VC_CM_array  = new double* [Nalpha*Nalpha];
	create_VC_col_maj_3N_pointer_array(VC_CM_array,
									   C_WP_unco_array,
									   C_WP_coup_array,
									   V_WP_unco_array,
									   V_WP_coup_array,
									   num_2N_unco_states,
									   num_2N_coup_states,
									   Np_WP,
									   pw_states,
									   run_parameters);
	
	size_t  dense_dim = Nalpha * Nq_WP * Np_WP;

	//// Temp code to print P123 rows
	//size_t row_idx = 5*Np_WP*Nq_WP + 1*Np_WP + 0;
	//for (size_t nnz_idx=0; nnz_idx<P123_sparse_dim; nnz_idx++){
	//	if (P123_sparse_row_array[nnz_idx]==row_idx){
	//		std::cout << P123_sparse_val_array[nnz_idx] << std::endl;
	//	}
	//}
	
	/* Test optimized routine for PVC columns */
	if (test_PVC_col_routine){
		printf("   - Testing PVC-column routine ... \n");
		PVC_col_calc_test(Nalpha,
						  Nq_WP,
						  Np_WP,
						  VC_CM_array,
						  P123_sparse_val_array,
						  P123_sparse_row_array,
						  P123_sparse_col_array_csc,
						  P123_sparse_dim);
		printf("     - Done \n");
	}
	/* Test optimized routine for CPVC columns */
	if (test_CPVC_col_routine){
		printf("   - Testing CPVC-column routine ... \n");
		CPVC_col_calc_test(Nalpha,
						   Nq_WP,
						   Np_WP,
						   CT_RM_array,
						   VC_CM_array,
						   P123_sparse_val_array,
						   P123_sparse_row_array,
						   P123_sparse_col_array_csc,
						   P123_sparse_dim);
		printf("     - Done \n");
	}
	
	auto timestamp_solve_start = std::chrono::system_clock::now();
	if (solve_dense==false){
		pade_method_solve(U_array,
						  U_BU_array,
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
						  P123_sparse_dim,
						  chn_os_indexing,
					      run_parameters,
						  file_identification);
	}
	else{
		printf("     - Solving Faddeev equation using a dense direct solver (WARNING: CAN TAKE LONG) ... \n");
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
	}

	auto timestamp_solve_end = std::chrono::system_clock::now();
	std::chrono::duration<double> time_solve = timestamp_solve_end - timestamp_solve_start;
	printf("     - Done. Time used: %.6f\n", time_solve.count());
}

/* Create array of pointers to C^T matrices for product (C^T)PVC in row-major format */
void create_CT_row_maj_3N_pointer_array(double** CT_RM_array,
										double*  C_WP_unco_array,
										double*  C_WP_coup_array,
										int  	 num_2N_unco_states,
										int  	 num_2N_coup_states,
										size_t   Np_WP,
										pw_3N_statespace pw_states,
										run_params run_parameters){
	
	size_t Nalpha		  = pw_states.Nalpha;
	int*   L_2N_array	  = pw_states.L_2N_array;
	int*   S_2N_array	  = pw_states.S_2N_array;
	int*   J_2N_array	  = pw_states.J_2N_array;
	int*   T_2N_array	  = pw_states.T_2N_array;
	int*   L_1N_array	  = pw_states.L_1N_array;
	int*   two_J_1N_array = pw_states.two_J_1N_array;
	int*   two_T_3N_array = pw_states.two_T_3N_array;

	/* This test will be reused several times */
	bool tensor_force_true = (run_parameters.tensor_force==true);
	
	/* Number of uncoupled and coupled 2N-channels */
	int num_unco_chns = num_2N_unco_states;
	int num_coup_chns = num_2N_coup_states;

	double* C_subarray  = NULL;
	double* CT_subarray = NULL;

	double* CT_unco_array = new double [Np_WP*Np_WP   * num_unco_chns];
	double* CT_coup_array = new double [Np_WP*Np_WP*4 * num_coup_chns];
	
	/* Copy and transpose all 2N-uncoupled C-arrays */
	for (int idx_chn_unco=0; idx_chn_unco<num_unco_chns; idx_chn_unco++){
		size_t idx_2N_mat_WP_unco = idx_chn_unco*Np_WP*Np_WP;
			
		C_subarray  = &C_WP_unco_array[idx_2N_mat_WP_unco];
		CT_subarray = &CT_unco_array  [idx_2N_mat_WP_unco];

		/* Copy content to avoid rewriting C-arrays */
		std::copy(C_subarray, C_subarray + Np_WP*Np_WP, CT_subarray);
			
		/* Transpose C to get C^T */
		simple_transpose_matrix_routine(CT_subarray, Np_WP);

		//double* tempprod = new double [Np_WP*Np_WP];
		//dot_MM(CT_subarray, C_subarray, tempprod, Np_WP, Np_WP, Np_WP);
		//for (int i=0; i<Np_WP; i++){
		//	for (int j=0; j<Np_WP; j++){
		//		if (i!=j and abs(tempprod[i*Np_WP+j])>1e-15){
		//			std::cout << "unco " << i << " " << j << " " << tempprod[i*Np_WP+j] << std::endl;
		//		}
		//	}
		//}
	}

	/* Copy and transpose all 2N-coupled C-arrays */
	for (int idx_chn_coup=0; idx_chn_coup<num_coup_chns; idx_chn_coup++){
		size_t idx_2N_mat_WP_coup = idx_chn_coup*4*Np_WP*Np_WP;

		C_subarray  = &C_WP_coup_array[idx_2N_mat_WP_coup];
		CT_subarray = &CT_coup_array  [idx_2N_mat_WP_coup];
		
		/* Copy content to avoid rewriting C-arrays */
		std::copy(C_subarray, C_subarray + 4*Np_WP*Np_WP, CT_subarray);
		
		/* Transpose C to get C^T */
		simple_transpose_matrix_routine(CT_subarray, 2*Np_WP);

		//double* tempprod = new double [4*Np_WP*Np_WP];
		//dot_MM(CT_subarray, C_subarray, tempprod, 2*Np_WP, 2*Np_WP, 2*Np_WP);
		//for (int i=0; i<2*Np_WP; i++){
		//	for (int j=0; j<2*Np_WP; j++){
		//		if (i!=j and abs(tempprod[i*2*Np_WP+j])>1e-15){
		//			std::cout << "coup " << i << " " << j << " " << tempprod[i*2*Np_WP+j] << std::endl;
		//		}
		//	}
		//}

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
		int two_T_3N_r = two_T_3N_array[idx_alpha_r];

		/* Column state */
		for (size_t idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
			int L_2N_c 	   = L_2N_array[idx_alpha_c];
			int S_2N_c 	   = S_2N_array[idx_alpha_c];
			int J_2N_c 	   = J_2N_array[idx_alpha_c];
			int T_2N_c 	   = T_2N_array[idx_alpha_c];
			int L_1N_c 	   = L_1N_array[idx_alpha_c];
			int two_J_1N_c = two_J_1N_array[idx_alpha_c];
			int two_T_3N_c = two_T_3N_array[idx_alpha_c];

			/* Check if possible channel through interaction */
			bool check_T = (T_2N_r==T_2N_c);
			bool check_J = (J_2N_r==J_2N_c);
			bool check_S = (S_2N_r==S_2N_c);
			bool check_L = ( (tensor_force_true && abs(L_2N_r-L_2N_c)<=2) || L_2N_r==L_2N_c);
			bool check_l = (L_1N_r==L_1N_c);
			bool check_j = (two_J_1N_r==two_J_1N_c);

			/* Check if possible channel through interaction */
			if (check_T && check_J && check_S && check_L && check_l && check_j){

				/* Detemine if this is a coupled channel.
				 * !!! With isospin symmetry-breaking we count 1S0 as a coupled matrix via T_3N-coupling !!! */
				bool coupled_matrix = false;
				bool state_1S0 = (S_2N_r==0 && J_2N_r==0 && L_2N_r==0);
				bool coupled_via_L_2N = (tensor_force_true && (L_2N_r!=L_2N_c || (L_2N_r==L_2N_c && L_2N_r!=J_2N_r && J_2N_r!=0)));
				bool coupled_via_T_3N = (state_1S0==true && run_parameters.isospin_breaking_1S0==true);
				if (coupled_via_L_2N && coupled_via_T_3N){
					raise_error("Warning! Code has not been written to handle isospin-breaking in coupled channels!");
				}
				if (coupled_via_L_2N || coupled_via_T_3N){ // This counts 3P0 as uncoupled; used in matrix structure
					coupled_matrix  = true;
				}

				/* find which VC-product corresponds to the current coupling */
				if (coupled_matrix){
					size_t idx_chn_coup       = (size_t) unique_2N_idx(L_2N_r, S_2N_r, J_2N_r, T_2N_r, coupled_matrix, run_parameters);
					size_t idx_2N_mat_WP_coup = idx_chn_coup*4*Np_WP*Np_WP;
					if ( (coupled_via_L_2N && L_2N_r<L_2N_c) || (coupled_via_T_3N && two_T_3N_r<two_T_3N_c) ){       // L_r=J_r-1, L_c=J_r+1 OR (for 1S0) two_T_3N_r==1/2, two_T_3N_c==3/2
						CT_subarray = &CT_coup_array[idx_2N_mat_WP_coup + 1*Np_WP*Np_WP];
					}
					else if ( (coupled_via_L_2N && L_2N_r>L_2N_c) || (coupled_via_T_3N && two_T_3N_r>two_T_3N_c) ){  // L_r=J_r+1, L_c=J_r-1 OR (for 1S0) two_T_3N_r==3/2, two_T_3N_c==1/2
						CT_subarray = &CT_coup_array[idx_2N_mat_WP_coup + 2*Np_WP*Np_WP];
					}
					else if ( (coupled_via_L_2N && L_2N_r<J_2N_c) || (coupled_via_T_3N && two_T_3N_r==1) ){ 		 // L_r=J_r-1, L_c=J_r-1 OR (for 1S0) two_T_3N_r==1/2, two_T_3N_c==1/2
						CT_subarray = &CT_coup_array[idx_2N_mat_WP_coup + 0*Np_WP*Np_WP];
					}
					else if ( (coupled_via_L_2N && L_2N_r>J_2N_c) || (coupled_via_T_3N && two_T_3N_r==3) ){  		 // L_r=J_r+1, L_c=J_r+1 OR (for 1S0) two_T_3N_r==3/2, two_T_3N_c==3/2
						CT_subarray = &CT_coup_array[idx_2N_mat_WP_coup + 3*Np_WP*Np_WP];
					}
					else{
						raise_error("Unknown 2N coupled-matrix requested in CT-array setup.");
					}
				}
				else{
					size_t idx_chn_unco       = (size_t) unique_2N_idx(L_2N_r, S_2N_r, J_2N_r, T_2N_r, coupled_matrix, run_parameters);
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
										int  	 num_2N_unco_states,
										int  	 num_2N_coup_states,
										size_t   Np_WP,
										pw_3N_statespace pw_states,
										run_params run_parameters){
	
	size_t Nalpha		  = pw_states.Nalpha;
	int*   L_2N_array	  = pw_states.L_2N_array;
	int*   S_2N_array	  = pw_states.S_2N_array;
	int*   J_2N_array	  = pw_states.J_2N_array;
	int*   T_2N_array	  = pw_states.T_2N_array;
	int*   L_1N_array	  = pw_states.L_1N_array;
	int*   two_J_1N_array = pw_states.two_J_1N_array;
	int*   two_T_3N_array = pw_states.two_T_3N_array;
	
	/* This test will be reused several times */
	bool tensor_force_true = (run_parameters.tensor_force==true);

	/* Number of uncoupled and coupled 2N-channels */
	int num_unco_chns = num_2N_unco_states;
	int num_coup_chns = num_2N_coup_states;

	double* V_subarray = NULL;
	double* C_subarray = NULL;

	double* VC_unco_array = new double [Np_WP*Np_WP   * num_unco_chns];
	double* VC_coup_array = new double [Np_WP*Np_WP*4 * num_coup_chns];

	double* VC_product = NULL;

	/* Calculate all 2N-uncoupled VC-products and convert to column-major format */
	for (int idx_chn_unco=0; idx_chn_unco<num_unco_chns; idx_chn_unco++){
		size_t idx_2N_mat_WP_unco = idx_chn_unco*Np_WP*Np_WP;
			
		V_subarray = &V_WP_unco_array[idx_2N_mat_WP_unco];
		C_subarray = &C_WP_unco_array[idx_2N_mat_WP_unco];
		VC_product = &VC_unco_array  [idx_2N_mat_WP_unco];

		/* Multiply V and C using BLAS */
		dot_MM(V_subarray, C_subarray, VC_product, Np_WP, Np_WP, Np_WP);

		/* Transpose VC-product to get column-major format */
		simple_transpose_matrix_routine(VC_product, Np_WP);
	}

	/* Calculate all 2N-coupled VC-products and convert to column-major format */
	for (int idx_chn_coup=0; idx_chn_coup<num_coup_chns; idx_chn_coup++){
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
		int two_T_3N_r = two_T_3N_array[idx_alpha_r];

		/* Column state */
		for (size_t idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
			int L_2N_c     = L_2N_array[idx_alpha_c];
			int S_2N_c     = S_2N_array[idx_alpha_c];
			int J_2N_c     = J_2N_array[idx_alpha_c];
			int T_2N_c     = T_2N_array[idx_alpha_c];
			int L_1N_c 	   = L_1N_array[idx_alpha_c];
			int two_J_1N_c = two_J_1N_array[idx_alpha_c];
			int two_T_3N_c = two_T_3N_array[idx_alpha_c];

			/* Check if possible channel through interaction */
			bool check_T = (T_2N_r==T_2N_c);
			bool check_J = (J_2N_r==J_2N_c);
			bool check_S = (S_2N_r==S_2N_c);
			bool check_L = ( (tensor_force_true && abs(L_2N_r-L_2N_c)<=2) || L_2N_r==L_2N_c);
			bool check_l = (L_1N_r==L_1N_c);
			bool check_j = (two_J_1N_r==two_J_1N_c);

			/* Check if possible channel through interaction */
			if (check_T && check_J && check_S && check_L && check_l && check_j){

				/* Detemine if this is a coupled channel.
				 * !!! With isospin symmetry-breaking we count 1S0 as a coupled matrix via T_3N-coupling !!! */
				bool coupled_matrix = false;
				bool state_1S0 = (S_2N_r==0 && J_2N_r==0 && L_2N_r==0);
				bool coupled_via_L_2N = (tensor_force_true && (L_2N_r!=L_2N_c || (L_2N_r==L_2N_c && L_2N_r!=J_2N_r && J_2N_r!=0)));
				bool coupled_via_T_3N = (state_1S0==true && run_parameters.isospin_breaking_1S0==true);
				if (coupled_via_L_2N && coupled_via_T_3N){
					raise_error("Warning! Code has not been written to handle isospin-breaking in coupled channels!");
				}
				if (coupled_via_L_2N || coupled_via_T_3N){ // This counts 3P0 as uncoupled; used in matrix structure
					coupled_matrix  = true;
				}

				/* find which VC-product corresponds to the current coupling */
				if (coupled_matrix){
					size_t idx_chn_coup       = (size_t) unique_2N_idx(L_2N_r, S_2N_r, J_2N_r, T_2N_r, coupled_matrix, run_parameters);
					size_t idx_2N_mat_WP_coup = idx_chn_coup*4*Np_WP*Np_WP;
					if ( (coupled_via_L_2N && L_2N_r<L_2N_c) || (coupled_via_T_3N && two_T_3N_r<two_T_3N_c) ){       // L_r=J_r-1, L_c=J_r+1 OR (for 1S0) two_T_3N_r==1/2, two_T_3N_c==3/2
						VC_product = &VC_coup_array[idx_2N_mat_WP_coup + 1*Np_WP*Np_WP];
					}
					else if ( (coupled_via_L_2N && L_2N_r>L_2N_c) || (coupled_via_T_3N && two_T_3N_r>two_T_3N_c) ){  // L_r=J_r+1, L_c=J_r-1 OR (for 1S0) two_T_3N_r==3/2, two_T_3N_c==1/2
						VC_product = &VC_coup_array[idx_2N_mat_WP_coup + 2*Np_WP*Np_WP];
					}
					else if ( (coupled_via_L_2N && L_2N_r<J_2N_c) || (coupled_via_T_3N && two_T_3N_r==1) ){ 		 // L_r=J_r-1, L_c=J_r-1 OR (for 1S0) two_T_3N_r==1/2, two_T_3N_c==1/2
						VC_product = &VC_coup_array[idx_2N_mat_WP_coup + 0*Np_WP*Np_WP];
					}
					else if ( (coupled_via_L_2N && L_2N_r>J_2N_c) || (coupled_via_T_3N && two_T_3N_r==3) ){  		 // L_r=J_r+1, L_c=J_r+1 OR (for 1S0) two_T_3N_r==3/2, two_T_3N_c==3/2
						VC_product = &VC_coup_array[idx_2N_mat_WP_coup + 3*Np_WP*Np_WP];
					}
					else{
						raise_error("Unknown 2N coupled-matrix requested in VC-array setup.");
					}
				}
				else{
					size_t idx_chn_unco       = (size_t) unique_2N_idx(L_2N_r, S_2N_r, J_2N_r, T_2N_r, coupled_matrix, run_parameters);
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
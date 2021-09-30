
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
									col_array[idx_CPVC] += prod;

									//int nnz_idx = row_to_nnz_array[idx_CPVC];
									//if (nnz_idx==-1){
									//	row_to_nnz_array[idx_CPVC] = num_nnz;
									//	nnz_to_row_array[num_nnz]  = idx_CPVC;
									//	nnz_idx = num_nnz;
									//	num_nnz += 1;
									//}
									//col_array[nnz_idx] += prod;
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


void calculate_PVC_col_gemm(double*  PVC_col_chunk,
							size_t	 idx_col_start,
							size_t	 idx_col_end,
							size_t   chunk_steplength,
							double*  P_col_chunk,
					        size_t   Nalpha,
							size_t   Nq_WP,
							size_t   Np_WP,
					        double** VC_CM_array,
					        double*  P123_val_array,
					        int*     P123_row_array,
					        size_t*  P123_col_array,
					        size_t   P123_dim){

	/* Dense dimension */
	size_t dense_dim = Nalpha*Nq_WP*Np_WP;

	/* Number of columns in chunk */
	size_t chunk_dim = idx_col_end - idx_col_start;

	/* Ensure PVC-chunk contains only zeroes */
	for (size_t idx=0; idx<dense_dim*chunk_steplength; idx++){
		PVC_col_chunk[idx] = 0;
	}

	double* VC_subarray = NULL;

	/* Loop through columns in steps of Np_WP*Nq_WP */
	size_t idx_alpha_c_start = idx_col_start / (Np_WP*Nq_WP);
	size_t idx_alpha_c_end   = idx_col_end   / (Np_WP*Nq_WP);

	for (size_t idx_alpha_c=idx_alpha_c_start; idx_alpha_c<idx_alpha_c_end; idx_alpha_c++){

		/* Loop through columns in steps of Np_WP */
		size_t idx_q_c_start = (idx_col_start - idx_alpha_c*Nq_WP*Np_WP) / Np_WP;
		size_t idx_q_c_end   = (idx_col_end   - idx_alpha_c*Nq_WP*Np_WP) / Np_WP;

		for (size_t idx_q_c=idx_q_c_start; idx_q_c<idx_q_c_end; idx_q_c++){

			/* Ensure P-chunk contains only zeroes */
			for (size_t idx=0; idx<Np_WP*dense_dim; idx++){
				P_col_chunk[idx] = 0;
			}

			/* Fill P-chunk of size Np_WP*dense_dim */
			for (size_t idx_p_c=0; idx_p_c<Np_WP; idx_p_c++){
				size_t idx_col = idx_alpha_c*Nq_WP*Np_WP + idx_q_c*Np_WP + idx_p_c;
				size_t idx_nnz_start = P123_col_array[idx_col];
				size_t idx_nnz_end   = P123_col_array[idx_col];
				for (size_t idx_nnz=idx_nnz_start; idx_nnz<idx_nnz_end; idx_nnz++){
					size_t idx_row = P123_row_array[idx_nnz];

					/* NOTE THAT P = P123 + P132 = 2*P123 */
					P_col_chunk[idx_row*Np_WP + idx_col] = 2*P123_val_array[idx_nnz];
				}
			}

			for (size_t idx_alpha_j=0; idx_alpha_j<Nalpha; idx_alpha_j++){
				VC_subarray = VC_CM_array[idx_alpha_c*Nalpha + idx_alpha_j];

				/* Only do inner-product if VC is not zero due to conservation laws */
				if (VC_subarray!=NULL){

					/* Multiply P-chunk with VC_subarray and fill corresponding part of PVC_col_chunk */
					double beta  = 1;
					double alpha = 1;
					MKL_INT M    = dense_dim;
					MKL_INT N    = Np_WP;
					MKL_INT K    = Np_WP;
					MKL_INT lda  = Np_WP;
					MKL_INT ldb  = Np_WP;
					MKL_INT ldc  = chunk_steplength;
					double* A    = P_col_chunk;
					double* B    = VC_subarray;
					double* C    = &PVC_col_chunk[idx_alpha_c*Nq_WP*Np_WP + idx_q_c*Np_WP - idx_col_start];
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);	// real multiplication
				}
			}
		}
	}
}
void calculate_CPVC_col_gemm(double*  CPVC_col_chunk,
							 size_t	  col_idx_start,
							 size_t	  col_idx_end,
							 size_t   chunk_steplength,
							 size_t   Nalpha,
							 size_t   Nq_WP,
							 size_t   Np_WP,
							 double** CT_RM_array,
							 double** VC_CM_array,
							 double*  P123_val_array,
							 int*     P123_row_array,
							 size_t*  P123_col_array,
							 size_t   P123_dim){
	
	/* Dense dimension */
	size_t dense_dim = Nalpha*Nq_WP*Np_WP;

	/* Number of columns in chunk */
	size_t chunk_dim = col_idx_end - col_idx_start;
	
	/* Generate PVC-column chunk */
	double* PVC_col_chunk = new double [dense_dim * chunk_steplength];

	/* Generate P-column chunk */
	double* P_col_chunk = new double [dense_dim * Np_WP];

	double* PVC_col = new double [dense_dim];
	for (size_t idx_col=col_idx_start; idx_col<col_idx_end; idx_col++){
		size_t idx_alpha_c = idx_col / (Np_WP*Nq_WP);
		size_t idx_q_c     = (idx_col % (Np_WP*Nq_WP)) /  Np_WP;
		size_t idx_p_c     = idx_col %  Np_WP;
		calculate_PVC_col(PVC_col,
						  idx_alpha_c, idx_p_c, idx_q_c,
						  Nalpha,      Nq_WP,   Np_WP,
						  VC_CM_array,
						  P123_val_array,
						  P123_row_array,
						  P123_col_array,
						  P123_dim);
		for (size_t idx=0; idx<dense_dim; idx++){
			PVC_col_chunk[idx*chunk_steplength + idx_col-col_idx_start] = PVC_col[idx];
		}
	}
	delete [] PVC_col;

	//calculate_PVC_col_gemm(PVC_col_chunk,
	//					   col_idx_start,
	//					   col_idx_end,
	//					   chunk_steplength,
	//					   P_col_chunk,
	//				       Nalpha,
	//					   Nq_WP,
	//					   Np_WP,
	//				       VC_CM_array,
	//				       P123_val_array,
	//				       P123_row_array,
	//				       P123_col_array,
	//				       P123_dim);

	double* CT_subarray = NULL;

	/* Loop over rows of CPVC by alpha and q */
	for (size_t idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
		for (size_t idx_q_r=0; idx_q_r<Nq_WP; idx_q_r++){
			/* Loop over alpha in inner-product C^T x PVC */
			for (size_t idx_alpha_i=0; idx_alpha_i<Nalpha; idx_alpha_i++){
				
				CT_subarray = CT_RM_array[idx_alpha_r*Nalpha + idx_alpha_i];
				
				/* Only do inner-product if C^T is not zero due to conservation laws */
				if (CT_subarray!=NULL){
					/* Multiply PVC-chunk with CT_subarray and fill corresponding part of CPVC_col_chunk */
					double beta  = 1;
					double alpha = 1;
					MKL_INT M    = Np_WP;
					MKL_INT N    = chunk_dim;
					MKL_INT K    = Np_WP;
					MKL_INT lda  = Np_WP;
					MKL_INT ldb  = chunk_steplength;
					MKL_INT ldc  = chunk_steplength;
					double* A    = CT_subarray;
					double* B    = &PVC_col_chunk[ (idx_alpha_i*Nq_WP*Np_WP + idx_q_r*Np_WP)*chunk_steplength];
					double* C    = &CPVC_col_chunk[(idx_alpha_r*Nq_WP*Np_WP + idx_q_r*Np_WP)*chunk_steplength];
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);	// real multiplication
				}
			}
		}
	}
	
	delete [] PVC_col_chunk;
	delete [] P_col_chunk;
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

				printf(" - U-matrix element for alpha'=%d, alpha=%d, q=%d: %.10e + %.10ei \n", idx_alpha_NDOS_row, idx_alpha_NDOS_col, idx_q_NDOS, U_array[idx_NDOS].real(), U_array[idx_NDOS].imag());
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
					   size_t    P123_sparse_dim,
					   run_params run_parameters,
					   std::string file_identification){

	/* Print Pade-approximant convergences */
	bool print_PA_convergences = true;
	/* Print Neumann terms */
	bool print_neumann_terms   = true;
	/* Store Neumann terms */
	bool store_neumann_terms   = true;
	/* Store An-matrices */
	bool store_An_arrays 	   = false;

	/* Use as many threads as possible in MKL-GEMM */
	mkl_set_num_threads(omp_get_max_threads());

	//MKL_NUM_THREADS = omp_get_max_threads();
	//printf("   - Will run with %d MKL-threads \n",mkl_get_max_threads()); fflush(stdout);

	/* Number of on-shell nucleon-deuteron channels (deuteron states can mix, hence ^2) */
	size_t num_on_shell_A_rows = num_deuteron_states * num_q_com;
	size_t num_on_shell_A_vals = num_deuteron_states * num_deuteron_states * num_q_com;
	
	/* Upper limit on polynomial approximation of Faddeev eq. */
	size_t NM_max = 14;
	size_t num_neumann_terms = 2*NM_max+1;

	/* Coefficients for calculating Pade approximant */
	cdouble* a_coeff_array = new cdouble [ num_neumann_terms * num_on_shell_A_vals];

	/* Dense dimension of 3N-channel */
	size_t dense_dim = Nalpha * Nq_WP * Np_WP;

	/* Allocate row-arrays for A*A^n, where A=(C^T)(P)(VC) */
	double* re_A_An_row_array 	   = new double [dense_dim * num_on_shell_A_rows];
	double* re_A_An_row_array_prev = new double [dense_dim * num_on_shell_A_rows];
	double* im_A_An_row_array 	   = new double [dense_dim * num_on_shell_A_rows];
	double* im_A_An_row_array_prev = new double [dense_dim * num_on_shell_A_rows];
	
	/* Set A_An-arrays to zero */
	for (size_t i=0; i<dense_dim*num_on_shell_A_rows; i++){
		re_A_An_row_array[i] 	  = 0;
		re_A_An_row_array_prev[i] = 0;
		im_A_An_row_array[i] 	  = 0;
		im_A_An_row_array_prev[i] = 0;
	}

	/* File-paths for storing A_An_row_array and on-shell neumann terms */
	std::string A_An_row_filename      = run_parameters.output_folder + "/An_rows" + file_identification + ".txt";
	std::string neumann_terms_filename = run_parameters.output_folder + "/neumann_terms" + file_identification + ".txt";


	/* Set initial values for A_Kn_row_array, where K^n=1 for n=0 */
	printf("   - Working on Pade approximant P[N,M] for N=%d, M=%d \n",0,0); fflush(stdout);
	printf("     - Calculating on-shell rows of A*K^n for n=%d. \n", 0); fflush(stdout);
	auto timestamp_start = std::chrono::system_clock::now();
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
	auto timestamp_end = std::chrono::system_clock::now();
	std::chrono::duration<double> time = timestamp_end - timestamp_start;
	printf("       - Time generating CPVC-rows:     %.6f \n", time.count()); fflush(stdout);
	printf("       - Done \n", time.count()); fflush(stdout);
	
	/* First Neumann-term */
	printf("     - Extracting on-shell Neumann-series terms a_n=A*K^n for n=%d. \n",0); fflush(stdout);
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
				cdouble a_coeff = re_A_An_row_array_prev[idx_row_NDOS*dense_dim + idx_col_NDOS];

				/* Store coefficient */
				size_t idx_NDOS = idx_d_row*num_deuteron_states*num_q_com + idx_d_col*num_q_com + idx_q_com;
				a_coeff_array[idx_NDOS*num_neumann_terms] = a_coeff;
				
				if (print_neumann_terms){
					printf("       - Neumann term %d for alpha'=%d, alpha=%d, q=%d: %.16e + %.16ei \n", 0, idx_alpha_NDOS_row, idx_alpha_NDOS_col, idx_q_NDOS, a_coeff.real(), a_coeff.imag());
					fflush(stdout);
				}
			}
		}
	}
	printf("       - Done \n"); fflush(stdout);

	if (store_An_arrays){
		printf("     - Storing matrix A*K^n for n=%d to output-folder. \n", 0); fflush(stdout);
		std::string array_seperator_text = "n = " + std::to_string(0);
		store_sep_complex_matrix(re_A_An_row_array_prev,
								 im_A_An_row_array_prev,
    	                         num_on_shell_A_rows,
						         dense_dim,
						         dense_dim,
						         A_An_row_filename,
						         true,
							     array_seperator_text);
		printf("       - Done \n");
	}
	if (store_neumann_terms){
		printf("     - Storing on-shell Neumann-series terms a_n=A*K^n for n=%d to output-folder. \n", 0); fflush(stdout);
		std::string array_seperator_text = "n = " + std::to_string(0);
		store_complex_matrix(a_coeff_array,
    	                     num_deuteron_states*num_deuteron_states*num_q_com,
			    			 num_neumann_terms,
			    			 num_neumann_terms,
			    			 neumann_terms_filename,
			    			 true,
							 array_seperator_text);
		printf("       - Done \n");
	}
	
	/* Arrays to store Pade-approximants (PA) for each on-shelle elements */
	cdouble* pade_approximants_array      = new cdouble [num_on_shell_A_vals * (NM_max+1)];
	size_t*  pade_approximants_idx_array  = new size_t  [num_on_shell_A_vals];
	bool*    pade_approximants_conv_array = new bool    [num_on_shell_A_vals];
	size_t	 num_converged_elements		  = 0;

	/* Define CPVC-chunks size */
	size_t num_bytes_per_chunk   = 0.5 * std::pow(1024,3);
    size_t num_bytes_in_CPVC_col = (sizeof(cdouble) * dense_dim);
    size_t num_cols_per_chunk    = num_bytes_per_chunk / num_bytes_in_CPVC_col;
	/* Ensure num_cols_per_chunk is divisible by Np_WP - important for gemm-use in CPVC-construction */
	num_cols_per_chunk			-= num_cols_per_chunk % Np_WP;
    size_t num_chunks            = dense_dim / num_cols_per_chunk + 1;
    size_t block_size 			 = num_cols_per_chunk;
	/* From test-script to program notation */
	size_t    max_num_cols_in_mem = block_size;
	size_t	  num_col_chunks	  = num_chunks;
	double*   CPVC_cols_array     = new double [dense_dim * max_num_cols_in_mem];
	
	/* Allocate row- and column-arrays for (C^T)(P)(VC) */
	int 	 num_threads			    = omp_get_max_threads();
	double*  omp_CPVC_col_array  		= new double [dense_dim * num_threads];
	int*     omp_CPVC_row_to_nnz_array  = new int    [dense_dim * num_threads];
	int*     omp_CPVC_nnz_to_row_array  = new int    [dense_dim * num_threads];

	for (size_t idx_NDOS=0; idx_NDOS<num_on_shell_A_vals; idx_NDOS++){
		pade_approximants_conv_array[idx_NDOS] = false;
	}

	/* Loop over number of Pade-terms we use */
	for (size_t NM=0; NM<NM_max+1; NM++){
		if (num_converged_elements==num_on_shell_A_vals){
			printf("   - Convergence reached for all on-shell elements! \n"); fflush(stdout);
			break;
		}

		if (NM!=0){
			printf("   - Working on Pade approximant P[N,M] for N=%d, M=%d \n",NM,NM); fflush(stdout);
		}
		
		size_t counter_array [100];
		for (size_t i=0; i<100; i++){
			counter_array[i] = 0;
		}

		/* Time-keeper array for parallel environment */
		double*  times_array = new double [3*num_threads];
		
		for (int n=2*NM-1; n<2*NM+1; n++){
			printf("     - Working on Neumann-terms for n=%d. \n", n); fflush(stdout);
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

			///* Reset time-keeper array */
			//auto timestamp_start = std::chrono::system_clock::now();
			//for (int i=0; i<3*num_threads; i++){
			//	times_array[i] = 0;
			//}

			/* Calculate all a-coefficients for calculated CPVC-column */
			double timestamp_resolvent_start = omp_get_wtime();
			printf("     - Multiplying in resolvent with An. \n"); fflush(stdout);
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
			
			printf("     - Calculating on-shell rows of A*K^n for n=%d. \n", n); fflush(stdout);
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

				/*void calculate_CPVC_col_gemm(double*  CPVC_col_chunk,
							 				 size_t	  col_idx_start,
							 				 size_t	  col_idx_end,
							 				 size_t   chunk_steplength,
							 				 int* 	  row_to_nnz_array, 
							 				 int* 	  nnz_to_row_array,
							 				 size_t&  num_nnz,
							 				 size_t   Nalpha,
							 				 size_t   Nq_WP,
							 				 size_t   Np_WP,
							 				 double** CT_RM_array,
							 				 double** VC_CM_array,
							 				 double*  P123_val_array,
							 				 int*     P123_row_array,
							 				 size_t*  P123_col_array,
							 				 size_t   P123_dim);*/
				calculate_CPVC_col_gemm(CPVC_cols_array,
										idx_col_start,
										idx_col_end,
										block_size,
										Nalpha, Nq_WP, Np_WP,
										CT_RM_array,
										VC_CM_array,
										P123_sparse_val_array,
										P123_sparse_row_array,
										P123_sparse_col_array,
										P123_sparse_dim);
										
				//#pragma omp parallel //num_threads(1)
				//{
				//size_t  thread_idx             = omp_get_thread_num();
				//double* CPVC_col_array  	   = &omp_CPVC_col_array	    [thread_idx*dense_dim];
				//int*    CPVC_row_to_nnz_array  = &omp_CPVC_row_to_nnz_array [thread_idx*dense_dim];
				//int*    CPVC_nnz_to_row_array  = &omp_CPVC_nnz_to_row_array [thread_idx*dense_dim];
				//#pragma omp for
				//for (size_t idx_col=idx_col_start; idx_col<idx_col_end; idx_col++){
				//	size_t idx_alpha_c = idx_col / (Np_WP*Nq_WP);
				//	size_t idx_q_c     = (idx_col % (Np_WP*Nq_WP)) /  Np_WP;
				//	size_t idx_p_c     = idx_col %  Np_WP;
				//
				//		/* Calculate CPVC-column */
				//		size_t CPVC_num_nnz = 0;
				//		size_t CPVC_col_idx = idx_col - idx_col_start;
				//		calculate_CPVC_col(&CPVC_cols_array[CPVC_col_idx*dense_dim],//CPVC_col_array,
				//						   CPVC_row_to_nnz_array,
				//						   CPVC_nnz_to_row_array,
				//						   CPVC_num_nnz,
				//						   idx_alpha_c, idx_p_c, idx_q_c,
				//						   Nalpha, Nq_WP, Np_WP,
				//						   CT_RM_array,
				//						   VC_CM_array,
				//						   P123_sparse_val_array,
				//						   P123_sparse_row_array,
				//						   P123_sparse_col_array,
				//						   P123_sparse_dim);
				//		counter_array[thread_idx] += CPVC_num_nnz;
				//}
				//}
				double timestamp_CPVC_chunk_end   = omp_get_wtime();
				time_CPVC_cols += timestamp_CPVC_chunk_end - timestamp_CPVC_chunk_start;

				//double beta  = 0;
				//double alpha = 1;
				//MKL_INT M   = num_on_shell_A_rows;
				//MKL_INT N   = cols_in_chunk;// max_num_cols_in_mem;
				//MKL_INT K   = dense_dim;
				//MKL_INT lda = dense_dim;
				//MKL_INT ldb = dense_dim;//max_num_cols_in_mem;
				//MKL_INT ldc = dense_dim;
				//double* re_A = &re_A_An_row_array_prev[0];
				//double* im_A = &im_A_An_row_array_prev[0];
				//double* B = &CPVC_cols_array[0];
				//double* re_C = &re_A_An_row_array[idx_col_start];
				//double* im_C = &im_A_An_row_array[idx_col_start];
				//
				//double timestamp_gemm_start = omp_get_wtime();
				//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, alpha, re_A, lda, B, ldb, beta, re_C, ldc);	// real multiplication
				//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, alpha, im_A, lda, B, ldb, beta, im_C, ldc);	// imag multiplication
				//double timestamp_gemm_end   = omp_get_wtime();
				//time_An_CPVC_multiply += timestamp_gemm_end - timestamp_gemm_start;

				double beta  = 0;
				double alpha = 1;
				MKL_INT M   = num_on_shell_A_rows;
				MKL_INT N   = cols_in_chunk;// max_num_cols_in_mem;
				MKL_INT K   = dense_dim;
				MKL_INT lda = dense_dim;
				MKL_INT ldb = max_num_cols_in_mem;//max_num_cols_in_mem;
				MKL_INT ldc = dense_dim;
				double* re_A = &re_A_An_row_array_prev[0];
				double* im_A = &im_A_An_row_array_prev[0];
				double* B = &CPVC_cols_array[0];
				double* re_C = &re_A_An_row_array[idx_col_start];
				double* im_C = &im_A_An_row_array[idx_col_start];
				
				double timestamp_gemm_start = omp_get_wtime();
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha, re_A, lda, B, ldb, beta, re_C, ldc);	// real multiplication
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha, im_A, lda, B, ldb, beta, im_C, ldc);	// imag multiplication
				double timestamp_gemm_end   = omp_get_wtime();
				time_An_CPVC_multiply += timestamp_gemm_end - timestamp_gemm_start;
				
				//double timestamp_4 = omp_get_wtime();
				//double time4 = timestamp_4 - timestamp_3;
				//times_array[3*0 + 2] += time4;
				
			}
			//auto timestamp_end = std::chrono::system_clock::now();
			//std::chrono::duration<double> time = timestamp_end - timestamp_start;
			double timestamp_neumann_end = omp_get_wtime();
			time_neumann = timestamp_neumann_end - timestamp_neumann_start;

			printf("       - Time multiplying An with G:    %.6f s \n", time_resolvent);
			printf("       - Time generating CPVC-cols:     %.6f s \n", time_CPVC_cols);
			printf("       - Time multiplying An with CPVC: %.6f s \n", time_An_CPVC_multiply);
			printf("       - Total time:                    %.6f s \n", time_neumann);
			printf("       - Done \n"); fflush(stdout);

			//double time_resetting_arrays = 0;
			//double time_CPVC_cols = 0;
			//double time_An_CPVC_multiply = 0;
			//for (int thread_idx=0; thread_idx<num_threads; thread_idx++){
			//	if (times_array[3*thread_idx]>time_resetting_arrays){
			//		time_resetting_arrays = times_array[3*thread_idx];
			//	}
			//	if (times_array[3*thread_idx+1]>time_CPVC_cols){
			//		time_CPVC_cols = times_array[3*thread_idx+1];
			//	}
			//	if (times_array[3*thread_idx+2]>time_An_CPVC_multiply){
			//		time_An_CPVC_multiply = times_array[3*thread_idx+2];
			//	}
			//}
			//printf("       - Time resetting arrays:         %.6f \n", time_resetting_arrays);
			//printf("       - Time generating CPVC-cols:     %.6f \n", time_CPVC_cols);
			//printf("       - Time multiplying with An-rows: %.6f \n", time_An_CPVC_multiply);
			//printf("       - Total time:                    %.6f \n", time.count());
			//printf("       - Done \n"); fflush(stdout);

			size_t nnz_counts = 0;
			for (size_t i=0; i<100; i++){
				nnz_counts += counter_array[i];
			}
			printf("NUMBER OF NNZ ELEMENTS IN A: %zu \n", nnz_counts);
			printf("NUMBER OF NNZ ELEMENTS IN P: %zu \n", P123_sparse_dim);

			/* Rewrite previous A_An with current A_An */
			for (size_t i=0; i<num_on_shell_A_rows*dense_dim; i++){
				re_A_An_row_array_prev[i] = re_A_An_row_array[i];
				im_A_An_row_array_prev[i] = im_A_An_row_array[i];
				//if (n==27){std::cout << A_An_row_array_prev[i] << std::endl;}
			}

			printf("     - Extracting on-shell Neumann-series terms a_n=A*K^n for n=%d. \n", n); fflush(stdout);
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
							printf("       - Neumann term %d for alpha'=%d, alpha=%d, q=%d: %.16e + %.16ei \n", n, idx_alpha_NDOS_row, idx_alpha_NDOS_col, idx_q_NDOS, a_coeff.real(), a_coeff.imag());
							//printf("%d\n", idx_row_NDOS*dense_dim + idx_col_NDOS);
						}
					}
				}
			}
			printf("       - Done \n");

			if (store_An_arrays){
			printf("     - Storing matrix A*K^n for n=%d to output-folder. \n", n); fflush(stdout);
			std::string array_seperator_text = "n = " + std::to_string(n);
		    store_sep_complex_matrix(re_A_An_row_array,
									 im_A_An_row_array,
                                     num_on_shell_A_rows,
		    				         dense_dim,
		    				         dense_dim,
		    				         A_An_row_filename,
		    				         false,
								     array_seperator_text);
			printf("       - Done \n");
			}
			if (store_neumann_terms){
				printf("     - Storing on-shell Neumann-series terms a_n=A*K^n for n=%d to output-folder. \n", n); fflush(stdout);
				std::string array_seperator_text = "n = " + std::to_string(n);
				store_complex_matrix(a_coeff_array,
            	                     num_deuteron_states*num_deuteron_states*num_q_com,
		    					     num_neumann_terms,
		    					     num_neumann_terms,
		    					     neumann_terms_filename,
		    					     false,
									 array_seperator_text);
				printf("       - Done \n");
			}
		}
		delete [] times_array;

		printf("     - Calculating Pade approximants PA[%d,%d]. \n", NM, NM); fflush(stdout);
		/* Calculate Pade approximants (PA) */
		for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
			for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
				for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
					///* Nucleon-deuteron on-shell (NDOS) indices
					// * (deuteron bound-state p-index is alwasy 0 due to eigenvalue ordering in SWP construction) */
					//size_t idx_alpha_NDOS_row = deuteron_idx_array[idx_d_row];
					//size_t idx_alpha_NDOS_col = deuteron_idx_array[idx_d_col];
					//size_t idx_q_NDOS 	   	  = q_com_idx_array[idx_q_com];

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

					/* Condition for convergence: 3 or more iterations past minimum PA, or last iteration NM=NM_max */
					if (NM-idx_best_PA>4 || NM==NM_max){
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
		printf("       - Done \n"); fflush(stdout);
	}

	printf("   - Extracting on-shell U-matrix elements \n"); fflush(stdout);
	/* Set on-shell U-matrix elements equal "best" PA */
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
				printf("     - U-matrix element for alpha'=%d, alpha=%d, q=%d: %.10e + %.10ei \n", idx_alpha_NDOS_row, idx_alpha_NDOS_col, idx_q_NDOS, U_array[idx_NDOS].real(), U_array[idx_NDOS].imag());
			}
		}
	}
	printf("     - Done \n"); fflush(stdout);

	printf("     - Storing on-shell Neumann-series terms a_n=A*K^n for all n to output-folder. \n"); fflush(stdout);
	store_complex_matrix(a_coeff_array,
                         num_deuteron_states*num_deuteron_states*num_q_com,
					     num_neumann_terms,
					     num_neumann_terms,
					     neumann_terms_filename,
					     true,
						 "Neumann terms");
	printf("       - Done \n");

	delete [] re_A_An_row_array;
	delete [] re_A_An_row_array_prev;
	delete [] im_A_An_row_array;
	delete [] im_A_An_row_array_prev;
	delete [] a_coeff_array;
	delete [] CPVC_cols_array;
	delete [] omp_CPVC_col_array;
	delete [] omp_CPVC_row_to_nnz_array;
	delete [] omp_CPVC_nnz_to_row_array;
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
							 int*      two_J_1N_array,
							 bool 	   tensor_force_true,
							 std::string file_identification,
					         run_params run_parameters){
	
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
									   tensor_force_true,
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
									   tensor_force_true,
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
	if (solve_dense==false){
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
						  P123_sparse_dim,
					      run_parameters,
						  file_identification);
	}
	else{
		printf("   - Solving Faddeev equation using a dense direct solver (WARNING: CAN TAKE LONG) ... \n");
		//auto timestamp_solve_start = std::chrono::system_clock::now();

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

		//auto timestamp_solve_end = std::chrono::system_clock::now();
		//std::chrono::duration<double> time_solve = timestamp_solve_end - timestamp_solve_start;
		//printf("   - Done. Time used: %.6f\n", time_solve.count());
	}

	auto timestamp_solve_end = std::chrono::system_clock::now();
	std::chrono::duration<double> time_solve = timestamp_solve_end - timestamp_solve_start;
	printf("   - Done. Time used: %.6f\n", time_solve.count());
}

/* Create array of pointers to C^T matrices for product (C^T)PVC in row-major format */
void create_CT_row_maj_3N_pointer_array(double** CT_RM_array,
										double*  C_WP_unco_array,
										double*  C_WP_coup_array,
										bool     tensor_force_true,
										size_t   Np_WP,
										int      J_2N_max,
										size_t   Nalpha,
										int*     L_2N_array,
										int*     S_2N_array,
										int*     J_2N_array,
										int*     T_2N_array,
							 			int*     L_1N_array, 
							 			int*     two_J_1N_array){
	
	/* Number of uncoupled and coupled 2N-channels */
	int num_unco_chns = 0;
	int num_coup_chns = 0;
	if (tensor_force_true){
		num_unco_chns = 2*(J_2N_max+1);
		num_coup_chns =    J_2N_max;
	}
	else{
		num_unco_chns = 4*J_2N_max + 2;
	}

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

		/* Column state */
		for (size_t idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
			int L_2N_c 	   = L_2N_array[idx_alpha_c];
			int S_2N_c 	   = S_2N_array[idx_alpha_c];
			int J_2N_c 	   = J_2N_array[idx_alpha_c];
			int T_2N_c 	   = T_2N_array[idx_alpha_c];
			int L_1N_c 	   = L_1N_array[idx_alpha_c];
			int two_J_1N_c = two_J_1N_array[idx_alpha_c];

			/* Check if possible channel through interaction */
			bool check_T = (T_2N_r==T_2N_c);
			bool check_J = (J_2N_r==J_2N_c);
			bool check_S = (S_2N_r==S_2N_c);
			bool check_L = ( (tensor_force_true && abs(L_2N_r-L_2N_c)<=2) || L_2N_r==L_2N_c);
			bool check_l = (L_1N_r==L_1N_c);
			bool check_j = (two_J_1N_r==two_J_1N_c);

			/* Check if possible channel through interaction */
			if (check_T && check_J && check_S && check_L && check_l && check_j){

				/* Detemine if this is a coupled channel */
				bool coupled_matrix = false;
				if ( tensor_force_true && (L_2N_r!=L_2N_c || (L_2N_r==L_2N_c & L_2N_r!=J_2N_r & J_2N_r!=0)) ){ // This counts 3P0 as uncoupled; used in matrix structure
					coupled_matrix  = true;
				}

				/* find which VC-product corresponds to the current coupling */
				if (coupled_matrix){
					size_t idx_chn_coup       = (size_t) unique_2N_idx(L_2N_r, S_2N_r, J_2N_r, T_2N_r, tensor_force_true, coupled_matrix);
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
					size_t idx_chn_unco       = (size_t) unique_2N_idx(L_2N_r, S_2N_r, J_2N_r, T_2N_r, tensor_force_true, coupled_matrix);
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
										bool     tensor_force_true,
										size_t   Np_WP,
										int      J_2N_max,
										size_t   Nalpha,
										int*     L_2N_array,
										int*     S_2N_array,
										int*     J_2N_array,
										int*     T_2N_array,
							 			int*     L_1N_array, 
							 			int*     two_J_1N_array){

	/* Number of uncoupled and coupled 2N-channels */
	int num_unco_chns = 0;
	int num_coup_chns = 0;
	if (tensor_force_true){
		num_unco_chns = 2*(J_2N_max+1);
		num_coup_chns =    J_2N_max;
	}
	else{
		num_unco_chns = 4*J_2N_max + 2;
	}

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

		/* Column state */
		for (size_t idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
			int L_2N_c     = L_2N_array[idx_alpha_c];
			int S_2N_c     = S_2N_array[idx_alpha_c];
			int J_2N_c     = J_2N_array[idx_alpha_c];
			int T_2N_c     = T_2N_array[idx_alpha_c];
			int L_1N_c 	   = L_1N_array[idx_alpha_c];
			int two_J_1N_c = two_J_1N_array[idx_alpha_c];

			/* Check if possible channel through interaction */
			bool check_T = (T_2N_r==T_2N_c);
			bool check_J = (J_2N_r==J_2N_c);
			bool check_S = (S_2N_r==S_2N_c);
			bool check_L = ( (tensor_force_true && abs(L_2N_r-L_2N_c)<=2) || L_2N_r==L_2N_c);
			bool check_l = (L_1N_r==L_1N_c);
			bool check_j = (two_J_1N_r==two_J_1N_c);

			/* Check if possible channel through interaction */
			if (check_T && check_J && check_S && check_L && check_l && check_j){

				/* Detemine if this is a coupled channel */
				bool coupled_matrix = false;
				if ( tensor_force_true && (L_2N_r!=L_2N_c || (L_2N_r==L_2N_c & L_2N_r!=J_2N_r & J_2N_r!=0)) ){ // This counts 3P0 as uncoupled; used in matrix structure
					coupled_matrix  = true;
				}

				/* find which VC-product corresponds to the current coupling */
				if (coupled_matrix){
					size_t idx_chn_coup       = (size_t) unique_2N_idx(L_2N_r, S_2N_r, J_2N_r, T_2N_r, tensor_force_true, coupled_matrix);
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
					size_t idx_chn_unco       = (size_t) unique_2N_idx(L_2N_r, S_2N_r, J_2N_r, T_2N_r, tensor_force_true, coupled_matrix);
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
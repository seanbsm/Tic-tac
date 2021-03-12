
#include "make_permutation_matrix.h"

/* Add this if 6j takes a while */
void precalculate_Wigner_6j_symbols(){
}

void calculate_Gtilde_array(double* Gtilde_array,
							double* Atilde_array,
							int  Nx, double* x_array,
							int  Nq, double* q_array,
							int  Np, double* p_array,
							int  Nalpha,
							int  L_max,
							int* L_2N_array,
							int* L_1N_array,
							int  two_J_3N){

	bool print_Gtilde_progress = false;

	long int fullsize = Np*Nq*Nx*Nalpha*(Nalpha - 1);
	long int counter = 0;
	int frac_n, frac_o=0;
	
	#pragma omp parallel
	{
		#pragma omp for

		for (MKL_INT64 p_idx=0; p_idx<Np; p_idx++){
			for (MKL_INT64 q_idx=0; q_idx<Nq; q_idx++){
				for (MKL_INT64 x_idx=0; x_idx<Nx; x_idx++){
					for (MKL_INT64 alpha=0; alpha<Nalpha; alpha++){
						for (MKL_INT64 alpha_p=0; alpha_p<Nalpha; alpha_p++){
							Gtilde_array[alpha*Nalpha*Np*Nq*Nx + alpha_p*Np*Nq*Nx + p_idx*Nq*Nx + q_idx*Nx + x_idx]
								= Gtilde_new (p_array[p_idx], q_array[q_idx], x_array[x_idx], alpha, alpha_p, Nalpha, L_max, L_2N_array, L_1N_array, Atilde_array, two_J_3N);

							counter += 1;
			
							if (print_Gtilde_progress){
								frac_n = (100*counter)/fullsize;
								if (frac_n>frac_o){std::cout << frac_n << "%" << std::endl; frac_o=frac_n;}
							}
						}
					}
				}
			}
		}
	}
}

void calculate_Gtilde_subarray(double* Gtilde_subarray,
							   double* Atilde_subarray,
							   int  Nx, double* x_array,
							   int  Nq, double* q_array,
							   int  Np, double* p_array,
							   int  L_2N, int L_2N_prime,
							   int  L_1N, int L_1N_prime,
							   int  two_J_3N){

	bool print_Gtilde_progress = false;

	long int fullsize = Np*Nq*Nx;
	long int counter = 0;
	int frac_n, frac_o=0;
	
	for (MKL_INT64 p_idx=0; p_idx<Np; p_idx++){
		for (MKL_INT64 q_idx=0; q_idx<Nq; q_idx++){
			for (MKL_INT64 x_idx=0; x_idx<Nx; x_idx++){
				Gtilde_subarray[p_idx*Nq*Nx + q_idx*Nx + x_idx]
					= Gtilde_subarray_new (p_array[p_idx],
										   q_array[q_idx],
										   x_array[x_idx],
										   L_2N, L_2N_prime,
										   L_1N, L_1N_prime,
										   Atilde_subarray,
										   two_J_3N);

				counter += 1;
			
				if (print_Gtilde_progress){
					frac_n = (100*counter)/fullsize;
					if (frac_n>frac_o){std::cout << frac_n << "%" << std::endl; frac_o=frac_n;}
				}
			}
		}
	}
}

void calculate_permutation_matrix_for_3N_channel(double** P123_val_dense_array,
												 double** P123_val_sparse_array,
												 int** P123_row_array,
												 int** P123_col_array,
												 int& P123_dim,
												 bool use_dense_format,
												 int  Nq, double* q_array, double* wq_array, int Np_per_WP, int Np_WP, double *p_array_WP_bounds,
												 int  Np, double* p_array, double* wp_array, int Nq_per_WP, int Nq_WP, double *q_array_WP_bounds,
												 int  Nx, double* x_array, double* wx_array,
												 int  Nalpha,
												 int* L_2N_array,
												 int* S_2N_array,
												 int* J_2N_array,
												 int* T_2N_array,
												 int* L_1N_array,
												 int* two_J_1N_array,
												 int  two_J_3N,
												 int  two_T_3N){
	
	/* Notation change */
	int two_J = two_J_3N;
	int two_T = two_T_3N;
	
	int Nx_Gtilde = Nx;
	int Jj_dim = Nalpha;
	int Np_3N = Np;
	int Nq_3N = Nq;

	int* L12_Jj    = L_2N_array;
	int* S12_Jj    = S_2N_array;
	int* J12_Jj    = J_2N_array;
	int* T12_Jj    = T_2N_array;
	int* l3_Jj     = L_1N_array;
	int* two_j3_Jj = two_J_1N_array;

	bool print_content = true;

	double* p_3N = p_array;
	double* q_3N = q_array;
	/* End of notation change */

	/* START OF OLD CODE SEGMENT WITH OLD VARIABLE-NOTATION */
	/* This code calculates the geometric function Gtilde_{alpha,alpha'}(p',q',x) as an array */

	// determine optimized Lmax: Lmax = max(get_L)+max(get_l)
	int max_L12 = 0;
	int max_l3 = 0;
	int max_J12 = 0;

	for (int alpha = 0; alpha <= Jj_dim - 1; alpha++){
		if (J12_Jj[alpha] > max_J12) max_J12 = J12_Jj[alpha];
		if (L12_Jj[alpha] > max_L12) max_L12 = L12_Jj[alpha];
		if (l3_Jj[alpha] > max_l3) max_l3 = l3_Jj[alpha];
	}

	// for F_local_matrix prestorage and F_interpolate
	int lmax = GSL_MAX_INT(max_l3, max_L12) + 3; // for C4 it is possible to couple l=lmax with THREE Y_{1}^{mu}
	int l_interpolate_max = l_interpolate_max = 2 * (lmax - 3) + 3;

	if (print_content){
		std::cout << "   - lmax = " << lmax << ", l_interpolate_max = " << l_interpolate_max << "\n";
	}

	int Lmax = max_L12 + max_l3;
	int two_jmax_SixJ = 2 * lmax; // do we need to prestore 6j??
	
	// for angular integration in Gtilde
	double x_Gtilde  [Nx_Gtilde];
	double wx_Gtilde [Nx_Gtilde];

	calc_gauss_points (x_Gtilde, wx_Gtilde, -1.0, 1.0, Nx_Gtilde);

	if (print_content){std::cout << "   - Nalpha    = " <<  Jj_dim << std::endl;}
	if (print_content){std::cout << "   - Np        = " <<  Np_3N << std::endl;}
	if (print_content){std::cout << "   - Nq        = " <<  Nq_3N << std::endl;}
	if (print_content){std::cout << "   - Nx_Gtilde = " <<  Nx_Gtilde << std::endl;}

	if (print_content){
		std::cout << "   - prestore SixJ...\n";
	}
	int SixJ_size = (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1);
	double *SixJ_array = new double[SixJ_size];
	if (SixJ_size < 0){
		raise_error("SixJ_array in make_permutation_matrix had negative size, likely an integer overflow. Check your dimensions.");
	}
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

	if (print_content){
		std::cout << "   - Started working on Atilde\n";
	}
	MKL_INT64 Atilde_N = Jj_dim * Jj_dim * (Lmax + 1);
	double *Atilde_store = new double[Atilde_N];
	for (MKL_INT64 i=0; i<Atilde_N; i++){
		Atilde_store[i] = 0.0;
	}
	for (int alpha = 0; alpha <= Jj_dim - 1; alpha++){
		for (int alphaprime = 0; alphaprime <= Jj_dim - 1; alphaprime++){
			for (int Ltotal = 0; Ltotal <= Lmax; Ltotal++){
				Atilde_store[alpha * Jj_dim * (Lmax + 1) + alphaprime * (Lmax + 1) + Ltotal] = Atilde (alpha, alphaprime, Ltotal, Jj_dim, L12_Jj, l3_Jj, J12_Jj, two_j3_Jj, S12_Jj, T12_Jj, two_J, two_T, SixJ_array, two_jmax_SixJ);
			}
		}
	}
	if (print_content){
		std::cout << "   - Finished working on Atilde\n";
	}

	/*if (print_content){
		std::cout << "   - Started working on Gtilde\n";
	}
	MKL_INT64 Gtilde_N = Jj_dim * Jj_dim * Np_3N * Nq_3N * Nx_Gtilde;
	double *Gtilde_store = new double[Gtilde_N];
	calculate_Gtilde_array(Gtilde_store,
						   Atilde_store,
						   Nx_Gtilde, x_Gtilde,
						   Nq_3N, q_3N,
						   Np_3N, p_3N,
						   Jj_dim,
						   Lmax,
						   L12_Jj,
						   l3_Jj,
						   two_J);
	if (print_content){
		std::cout << "   - Finished working on Gtilde\n";
	}*/

	/* END OF OLD CODE SEGMENT WITH OLD VARIABLE-NOTATION */

	/* Precalculate overlapping bins and where p_bar and q_bar are non-zero */
	bool* pq_WP_overlap_array = new bool [Nq_WP*Nq_WP*Np_WP*Np_WP];

	double qp_WP_bound_lower = 0; double q_WP_bound_lower = 0;
	double qp_WP_bound_upper = 0; double q_WP_bound_upper = 0;

	double pp_WP_bound_lower = 0; double p_WP_bound_lower = 0;
	double pp_WP_bound_upper = 0; double p_WP_bound_upper = 0;

	int counter = 0;
	for (int qp_idx_WP = 0; qp_idx_WP < Nq_WP; qp_idx_WP++){
		qp_WP_bound_lower = q_array_WP_bounds[qp_idx_WP];
		qp_WP_bound_upper = q_array_WP_bounds[qp_idx_WP+1];
		for (int pp_idx_WP = 0; pp_idx_WP < Np_WP; pp_idx_WP++){
			pp_WP_bound_lower = p_array_WP_bounds[pp_idx_WP];
			pp_WP_bound_upper = p_array_WP_bounds[pp_idx_WP+1];
			
			double p_bar_max = 0;
			double p_bar_min = 0;
			double q_bar_max = 0;
			double q_bar_min = 0;

			p_bar_max = pi1_tilde(pp_WP_bound_upper, qp_WP_bound_upper, +1); // x=+1
			p_bar_min = pi1_tilde(pp_WP_bound_lower, qp_WP_bound_lower, -1); // x=-1
			q_bar_max = pi2_tilde(pp_WP_bound_upper, qp_WP_bound_upper, -1); // x=-1
			q_bar_min = pi2_tilde(pp_WP_bound_lower, qp_WP_bound_lower, +1); // x=+1

			int pp_idx_lower = Np_per_WP* pp_idx_WP;
    		int pp_idx_upper = Np_per_WP*(pp_idx_WP+1);

    		int qp_idx_lower = Nq_per_WP* qp_idx_WP;
    		int qp_idx_upper = Nq_per_WP*(qp_idx_WP+1);

			for (int pp_idx=pp_idx_lower; pp_idx<pp_idx_upper; pp_idx++){
				double pp = p_array[pp_idx];
				if (pp<pp_WP_bound_lower || pp_WP_bound_upper<pp){
					raise_error("Found bin-quadrature p-momentum outside of bin-boundaries");
				}
        		for (int qp_idx=qp_idx_lower; qp_idx<qp_idx_upper; qp_idx++){
					double qp = q_array[qp_idx];
					if (qp<qp_WP_bound_lower || qp_WP_bound_upper<qp){
						raise_error("Found bin-quadrature q-momentum outside of bin-boundaries");
					}
					for (int x_idx=0; x_idx<Nx; x_idx++){
                		double x = x_array[x_idx];
						double p_bar = pi1_tilde(pp, qp, x);
						double q_bar = pi2_tilde(pp, qp, x);
						if (p_bar>p_bar_max){
							p_bar_max = p_bar;
						}
						if (p_bar<p_bar_min){
							p_bar_min = p_bar;
						}
						if (q_bar>q_bar_max){
							q_bar_max = q_bar;
						}
						if (q_bar<q_bar_min){
							q_bar_min = q_bar;
						}
					}
				}
			}

			if (p_bar_min>p_bar_max || q_bar_min>q_bar_max){
				raise_error("Wrongful p_bar or q_bar calculation in make_permutation_matrix.cpp");
			}

			for (int q_idx_WP = 0; q_idx_WP < Nq_WP; q_idx_WP++){
				q_WP_bound_lower = q_array_WP_bounds[q_idx_WP];
				q_WP_bound_upper = q_array_WP_bounds[q_idx_WP+1];
				for (int p_idx_WP = 0; p_idx_WP < Np_WP; p_idx_WP++){
					p_WP_bound_lower = p_array_WP_bounds[p_idx_WP];
					p_WP_bound_upper = p_array_WP_bounds[p_idx_WP+1];

					bool WP_overlap = true;

            		/* Check if p_bar can have values in current bin */
            		if ( p_bar_max<p_WP_bound_lower || p_WP_bound_upper<p_bar_min ){
						WP_overlap = false;
            		}
					/* Check if q_bar can have values in current bin */
            		if ( q_bar_max<q_WP_bound_lower || q_WP_bound_upper<q_bar_min ){
						WP_overlap = false;
            		}

					/* Unique index for current combination of WPs */
					int pq_WP_idx = qp_idx_WP*Np_WP*Nq_WP*Np_WP
							  	  + pp_idx_WP*Nq_WP*Np_WP
							  	  +  q_idx_WP*Np_WP
							  	  +  p_idx_WP;

					pq_WP_overlap_array[pq_WP_idx] = WP_overlap;

					if (WP_overlap==false){
						counter += 1;
					}
				}
			}
		}
	}
	printf("   - %.3f%% of P123-matrix violates momentum-conservation \n", 100.*counter/(Nq_WP*Nq_WP*Np_WP*Np_WP));

	/* Preallocate array if we use dense format. Otherwise (i.e. sparse) start with
	 * some reasonable guess (usually less than a percent), and expand if required. */
	long int P123_dense_dim    = Np_WP * Nq_WP * Nalpha;
	long int P123_dense_dim_sq = P123_dense_dim * P123_dense_dim;
	int sparse_step_length     = 0;

	if (use_dense_format){
		*P123_val_dense_array = new double [P123_dense_dim_sq];
		P123_dim        = P123_dense_dim;
	}
	else{
		if (P123_dense_dim_sq>1e6){
			sparse_step_length = P123_dense_dim;
		}
		else{
			sparse_step_length = P123_dense_dim_sq;
		}

		/* The sparse dimension is determined by counting, see the loops below */
		P123_dim = 0;

		*P123_val_sparse_array = new double [sparse_step_length];
		*P123_row_array 	   = new int    [sparse_step_length];
		*P123_col_array 	   = new int    [sparse_step_length];
	}
	
	int current_array_dim = sparse_step_length;

	int Gtilde_subarray_size = Np_per_WP * Nq_per_WP * Nx_Gtilde;

	double *P123_col_val_array = new double[P123_dense_dim];
	double  P123_col_tot_time  = 0;
	double  P123_av_row_time   = 0;
	int 	P123_rows_left	   = 0;
	double 	P123_est_time_left = 0;

	//int row_step_length = P123_dense_dim/100;
	//int row_multiplier = 1;
	
	/* <X_i'j'^alpha'| - loops (rows of P123) */
	for (int alphap_idx = 0; alphap_idx < Nalpha; alphap_idx++){
		for (int qp_idx_WP = 0; qp_idx_WP < Nq_WP; qp_idx_WP++){
			for (int pp_idx_WP = 0; pp_idx_WP < Np_WP; pp_idx_WP++){
				int P123_row_idx = alphap_idx*Nq_WP*Np_WP + qp_idx_WP*Np_WP +  pp_idx_WP;

				//if (P123_row_idx!=row_step_length*row_multiplier){
				//	continue;
				//}
				//else{
				//	row_multiplier += 1;
				//}

				/* Progress printout */
				if (print_content){
					P123_av_row_time = P123_col_tot_time/P123_row_idx;
					P123_rows_left   = P123_dense_dim - P123_row_idx;
					P123_est_time_left = P123_rows_left*P123_av_row_time / 3600.;
					printf("\r   - Working on row %d of %d. Av. time per row: %.3f s. Est. completion time: %.1f h", P123_row_idx, P123_dense_dim, P123_av_row_time, P123_est_time_left);
					fflush(stdout);
				}

				/* |X_ij^alpha> - loops (columns of P123) */
				auto timestamp_P123_col_start = chrono::system_clock::now();
				#pragma omp parallel
				{
					#pragma omp for
					for (int alpha_idx = 0; alpha_idx < Nalpha; alpha_idx++){

						int L_2N = L_2N_array[alphap_idx];
						int L_1N = L_1N_array[alphap_idx];

						int L_2N_prime = L_2N_array[alpha_idx];
						int L_1N_prime = L_1N_array[alpha_idx];

						double Gtilde_subarray [Gtilde_subarray_size];
						calculate_Gtilde_subarray(Gtilde_subarray,
												  &Atilde_store[alphap_idx*Nalpha*(Lmax+1) + alpha_idx*(Lmax+1)],
												  Nx, x_array,
												  Nq_per_WP, &q_array[qp_idx_WP*Nq_per_WP],
												  Np_per_WP, &p_array[pp_idx_WP*Np_per_WP],
												  L_2N, L_2N_prime,
												  L_1N, L_1N_prime,
												  two_J_3N);

						//printf("%d %d %d %d \n", alphap_idx, qp_idx_WP, pp_idx_WP, alpha_idx);
						//for (MKL_INT64 p_idx=pp_idx_WP*Np_per_WP; p_idx<(pp_idx_WP+1)*Np_per_WP; p_idx++){
						//	for (MKL_INT64 q_idx=qp_idx_WP*Nq_per_WP; q_idx<(qp_idx_WP+1)*Nq_per_WP; q_idx++){
						//		for (MKL_INT64 x_idx=0; x_idx<Nx; x_idx++){
						//			double G_sub  = Gtilde_subarray[(p_idx-pp_idx_WP*Np_per_WP)*Nq_per_WP*Nx + (q_idx-qp_idx_WP*Nq_per_WP)*Nx + x_idx];
						//			double G_full = Gtilde_store[alphap_idx*Nalpha*Np*Nq*Nx + alpha_idx*Np*Nq*Nx + p_idx*Nq*Nx + q_idx*Nx + x_idx];
						//			printf("%d %d %d %.16f %.16f \n", p_idx, q_idx, x_idx, G_sub, G_full);
						//			if (abs(G_full-G_sub)>1e-15){
						//				raise_error("fuck");
						//			}
						//		}
						//	}
						//}

						for (int q_idx_WP = 0; q_idx_WP < Nq_WP; q_idx_WP++){
							for (int p_idx_WP = 0; p_idx_WP < Np_WP; p_idx_WP++){

								int P123_col_idx = alpha_idx*Nq_WP*Np_WP + q_idx_WP*Np_WP + p_idx_WP;

								/* Unique index for current combination of WPs */
								int pq_WP_idx = qp_idx_WP*Np_WP*Nq_WP*Np_WP
										      + pp_idx_WP*Nq_WP*Np_WP
										      +  q_idx_WP*Np_WP
										      +  p_idx_WP;
								bool WP_overlap = pq_WP_overlap_array[pq_WP_idx];

								double P123_element = 0;
								/* Only calculate P123 if there is WP bin-overlap in Heaviside functions */
								if (WP_overlap){
									P123_element = calculate_P123_element_in_WP_basis (  alpha_idx,  p_idx_WP,  q_idx_WP, 
																					    alphap_idx, pp_idx_WP, qp_idx_WP, 
																					    Np_per_WP, p_array, wp_array,
																					    Nq_per_WP, q_array, wq_array,
																					    Nx,    x_array, wx_array,
																					    Np_WP, p_array_WP_bounds,
																					    Nq_WP, q_array_WP_bounds,
																					    Nalpha,
																					    Gtilde_subarray);
																					    //Gtilde_store );
								}

								if (use_dense_format){
									int P123_mat_idx = (int) P123_row_idx*P123_dense_dim + P123_col_idx;
									(*P123_val_dense_array)[P123_mat_idx] = P123_element;
								}
								else{
									P123_col_val_array[P123_col_idx] = P123_element;
								}
							}
						}
					}
				}

				/* Loop through column-element and append non-zero elements */
				for (int P123_col_idx = 0; P123_col_idx < P123_dense_dim; P123_col_idx++){
					double P123_element = P123_col_val_array[P123_col_idx];
					if (P123_element!=0){
						/* Append to sparse value and index arrays */
						(*P123_val_sparse_array)[P123_dim] = P123_element;
						(*P123_row_array)[P123_dim]        = P123_row_idx;
						(*P123_col_array)[P123_dim]        = P123_col_idx;

						/* Increment sparse dimension (num of non-zero elements) */
						P123_dim += 1;

						/* If the dimension goes over the array dimension we increase the array size
						 * via a copy-paste-type routine, and increment the current array dimension */
						if ( P123_dim>=current_array_dim ){  // This should occur a small amount of the time
							increase_sparse_array_size(P123_val_sparse_array, current_array_dim, sparse_step_length);
							increase_sparse_array_size(P123_row_array,        current_array_dim, sparse_step_length);
							increase_sparse_array_size(P123_col_array,        current_array_dim, sparse_step_length);

							/* Increment sparse-array dimension */
							current_array_dim += sparse_step_length;
						}
					}
				}

				auto timestamp_P123_col_end = chrono::system_clock::now();
				chrono::duration<double> time_P123_col = timestamp_P123_col_end - timestamp_P123_col_start;
				P123_col_tot_time += time_P123_col.count();
			}
		}
	}

	if (print_content){
		std::cout << "\n";
	}

	/* Contract arrays to minimal size (number of non-zero elements) */
	if (use_dense_format==false){
		reduce_sparse_array_size(P123_val_sparse_array, current_array_dim, P123_dim);
		reduce_sparse_array_size(P123_row_array,        current_array_dim, P123_dim);
		reduce_sparse_array_size(P123_col_array,        current_array_dim, P123_dim);
	}

	/* Delete all temporary arrays */
	//delete [] Gtilde_store;
	//delete [] Gtilde_subarray;
	delete [] pq_WP_overlap_array;
	delete [] Atilde_store;
	delete [] SixJ_array;
}

void calculate_permutation_matrices_for_all_3N_channels(double** P123_sparse_val_array,
														int**    P123_sparse_row_array,
														int**    P123_sparse_col_array,
														int&     P123_sparse_dim_array,
														int  Nq, double* q_array, double* wq_array, int Np_per_WP, int Np_WP, double *p_array_WP_bounds,
														int  Np, double* p_array, double* wp_array, int Nq_per_WP, int Nq_WP, double *q_array_WP_bounds,
														int  Nx, double* x_array, double* wx_array,
														int  Nalpha,
														int* L_2N_array,
														int* S_2N_array,
														int* J_2N_array,
														int* T_2N_array,
														int* L_1N_array,
														int* two_J_1N_array,
														int  two_J_3N,
														int  two_T_3N){
	
	bool print_content = true;

	/* Either use a dense-block or allocate space dynamically
	 * as non-zero elements are found. The 2nd approach is a bit slower,
	 * but necessary due to the typical size of the dense-block being unhandleable. */
	bool use_dense_format = false;
	
	/* Dense and sparse array dimensions */
	int     P123_array_dense_dim  = Nalpha * Np_WP * Nq_WP;
	int     P123_array_dim        = 0;

	/* Temporary dense matrix to be filled in subroutine if use_dense_format=true */
	double* P123_val_dense_array = NULL;

	if (print_content){
		printf(" - Begin calculating 3N-channel permutation matrix \n");
	}
	calculate_permutation_matrix_for_3N_channel(&P123_val_dense_array, P123_sparse_val_array, P123_sparse_row_array, P123_sparse_col_array, P123_array_dim, use_dense_format,
							  					Nq_WP*Nq_per_WP, q_array, wq_array, Np_per_WP, Np_WP, p_array_WP_bounds,
							  					Np_WP*Np_per_WP, p_array, wp_array, Nq_per_WP, Nq_WP, q_array_WP_bounds,
							  					Nx, x_array, wx_array,
							  					Nalpha,
							  					L_2N_array,
							  					S_2N_array,
							  					J_2N_array,
							  					T_2N_array,
							  					L_1N_array,
							  					two_J_1N_array,
												two_J_3N,
												two_T_3N);

	if (print_content){
		printf("   - Calculation finished. \n");
	}

	/* Convert dense-storage to sparse COO-format */
	if (use_dense_format){

		///* Old code snippet to check max matrix element, handy at times */ 
		//double max_element_dense = 0;
		//for (int idx=0; idx<P123_array_dim*P123_array_dim; idx++){
		//    if (abs(P123_val_dense_array[idx])>max_element_dense){
		//        max_element_dense = abs(P123_val_dense_array[idx]);
		//    } 
		//}
		//std::cout << "Max element dense: " << max_element_dense << std::endl;

		if (print_content){
			printf("   - Used dense format for calculation. Converting to sparse format \n");
		}
			
		square_dense_to_sparse_COO_format_converter(P123_array_dense_dim,
													P123_val_dense_array,
													P123_sparse_val_array,
													P123_sparse_row_array,
													P123_sparse_col_array,
													P123_sparse_dim_array);
	}
	else{
		P123_sparse_dim_array = P123_array_dim;
	}
		
	if (print_content){
		double P_123_subarray_density = (double) P123_sparse_dim_array / (P123_array_dense_dim*P123_array_dense_dim);

		printf(" - Dense dimension:   %dx%d \n", P123_array_dense_dim, P123_array_dense_dim);
		printf(" - Non-zero elements: %d \n",    P123_sparse_dim_array);
		printf(" - Density:           %.4f %% \n",  100*P_123_subarray_density);
	}
		
	delete [] P123_val_dense_array;
}
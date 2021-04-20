
#include "make_permutation_matrix.h"

/* Copied from: https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes */
template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  std::stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

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

void calculate_Gtilde_subarray_cart(double* Gtilde_subarray,
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

void calculate_Gtilde_subarray_polar(double* Gtilde_subarray,
								     double* Atilde_subarray,
								     int  Nx, double* x_array,
								     int  Nphi,
								     double* sin_phi_subarray,
								     double* cos_phi_subarray,
								     int  L_2N, int L_2N_prime,
								     int  L_1N, int L_1N_prime,
								     int  two_J_3N){

	bool print_Gtilde_progress = false;

	long int fullsize = Nphi*Nx;
	long int counter = 0;
	int frac_n, frac_o=0;
	
	for (MKL_INT64 phi_idx=0; phi_idx<Nphi; phi_idx++){
		for (MKL_INT64 x_idx=0; x_idx<Nx; x_idx++){
			Gtilde_subarray[phi_idx*Nx+ x_idx]
				= Gtilde_subarray_new (sin_phi_subarray[phi_idx],
									   cos_phi_subarray[phi_idx],
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

std::string generate_subarray_file_name(int two_J_3N, int two_T_3N, int P_3N,
										int Np_WP, int Nq_WP,
										int J_2N_max,
										int thread_idx,
										int current_TFC){

	std::string filename = "P123_subsparse_JTP_"
						 + to_string(two_J_3N) + "_" + to_string(two_T_3N) + "_" + to_string(P_3N)
						 + "_Np_" + to_string(Np_WP) + "_Nq_" + to_string(Nq_WP)
						 + "_J2max_" + to_string(J_2N_max) + "_TFC_" + to_string(thread_idx)
						 + "_" + to_string(current_TFC) + ".h5";

	return filename;
}

void calculate_permutation_matrix_for_3N_channel(double** P123_val_dense_array,
												 double** P123_val_sparse_array,
												 int**    P123_row_array,
												 int**    P123_col_array,
												 size_t&  P123_dim,
												 bool     use_dense_format,
												 bool     production_run,
												 int      Np_WP, double *p_array_WP_bounds,
												 int      Nq_WP, double *q_array_WP_bounds,
												 int      Nx, double* x_array, double* wx_array,
												 int      Nphi,
												 int      J_2N_max,
												 int      Nalpha,
												 int*     L_2N_array,
												 int*     S_2N_array,
												 int*     J_2N_array,
												 int*     T_2N_array,
												 int*     L_1N_array,
												 int*     two_J_1N_array,
												 int      two_J_3N,
												 int      two_T_3N,
												 int      P_3N){
	
	bool print_content = true;

	/* Notation change */
	int two_J = two_J_3N;
	int two_T = two_T_3N;
	
	int Nx_Gtilde = Nx;
	int Jj_dim = Nalpha;

	int* L12_Jj    = L_2N_array;
	int* S12_Jj    = S_2N_array;
	int* J12_Jj    = J_2N_array;
	int* T12_Jj    = T_2N_array;
	int* l3_Jj     = L_1N_array;
	int* two_j3_Jj = two_J_1N_array;

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
	if (print_content){std::cout << "   - Np_WP     = " <<  Np_WP << std::endl;}
	if (print_content){std::cout << "   - Nq_WP     = " <<  Nq_WP << std::endl;}
	if (print_content){std::cout << "   - Nphi      = " <<  Nphi << std::endl;}
	if (print_content){std::cout << "   - Nx        = " <<  Nx_Gtilde << std::endl;}

	if (print_content){
		std::cout << "   - prestore SixJ...\n";
	}
	int SixJ_size = (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1);
	double *SixJ_array = new double [SixJ_size];
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

	/* END OF OLD CODE SEGMENT WITH OLD VARIABLE-NOTATION */

	/* Precalculate overlapping bins and where p_bar and q_bar are non-zero */
	bool*    pq_WP_overlap_array = new bool   [Nq_WP*Nq_WP*Np_WP*Np_WP];
	double*   phi_array      	 = new double [Nq_WP*Np_WP*Nphi];
	double*  wphi_array      	 = new double [Nq_WP*Np_WP*Nphi];
	
	printf("   - Precalculating momentum conservations \n");
	fflush(stdout);
	#pragma omp parallel
	{
		#pragma omp for
	for (int qp_idx_WP=0; qp_idx_WP<Nq_WP; qp_idx_WP++){
		double qp_l = q_array_WP_bounds[qp_idx_WP];
		double qp_u = q_array_WP_bounds[qp_idx_WP+1];
		for (int pp_idx_WP=0; pp_idx_WP<Np_WP; pp_idx_WP++){
			double pp_l = p_array_WP_bounds[pp_idx_WP];
			double pp_u = p_array_WP_bounds[pp_idx_WP+1];
	
			double phi_lower = atan(pp_l/qp_u);
			double phi_upper = atan(pp_u/qp_l);
			
			/* Create phi-mesh */
			calc_gauss_points ( &phi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi],
							   &wphi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi],
							   phi_lower, phi_upper,
							   Nphi);
							   
			/* Verify that on-shell elements can exist for given phi-boundaries */
			for (int q_idx_WP=0; q_idx_WP<Nq_WP; q_idx_WP++){
				double q_l = q_array_WP_bounds[q_idx_WP];
				double q_u = q_array_WP_bounds[q_idx_WP+1];
				for (int p_idx_WP=0; p_idx_WP<Np_WP; p_idx_WP++){
					double p_l = p_array_WP_bounds[p_idx_WP];
					double p_u = p_array_WP_bounds[p_idx_WP+1];
	
					bool WP_overlap = false;
					/* Ensure possible phi boundaries */
					int hit_counter = 0;
					if (phi_lower<phi_upper){
						/* Search for on-shell elements */
						for (int phi_idx=0; phi_idx<Nphi; phi_idx++){
							double phi = phi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi + phi_idx];
							double sin_phi = std::sin(phi);
							double cos_phi = std::cos(phi);
	
							double kmin = std::max(pp_l/sin_phi, qp_l/cos_phi);
							double kmax = std::min(pp_u/sin_phi, qp_u/cos_phi);
	
							/* Ensure possible k-boundaries */
							if (kmin<kmax){
								for (int x_idx=0; x_idx<Nx; x_idx++){
									double x = x_array[x_idx];
	
									double zeta_1 = pi1_tilde(sin_phi, cos_phi, x);
									double zeta_2 = pi2_tilde(sin_phi, cos_phi, x);
	
									double Heaviside_lower = std::max(p_l/zeta_1, q_l/zeta_2);
									double Heaviside_upper = std::min(p_u/zeta_1, q_u/zeta_2);
	
									/* Ensure overlapping Heaviside boundaries */
									if (Heaviside_lower<Heaviside_upper){
										double kpmin = std::max(kmin, Heaviside_lower);
										double kpmax = std::min(kmax, Heaviside_upper);
	
										/* Skip momentum-violating integral boundaries */
										if (kpmin<kpmax){
											WP_overlap = true;
											
											//if (std::abs(kpmax-kpmin)>1e-15){
											//	hit_counter += 1;
											//}
											break;
										}
									}
								}
							}
							if (WP_overlap){
								break;
							}
						}
					}
	
					/* Unique index for current combination of WPs */
					size_t pq_WP_idx = qp_idx_WP*Np_WP*Nq_WP*Np_WP
							  	  	 + pp_idx_WP*Nq_WP*Np_WP
							  	  	 +  q_idx_WP*Np_WP
							  	  	 +  p_idx_WP;
	
					pq_WP_overlap_array[pq_WP_idx] = WP_overlap;
				}
			}
		}
	}
	}
	//double sparsity = 0.99881;
	//for (int i=0; i<Nq_WP*Nq_WP*Np_WP*Np_WP; i++){
	//	double prob_nnz = (double) rand() / RAND_MAX;
    //    if (prob_nnz>sparsity){
	//		pq_WP_overlap_array[i] = true;
	//	}
	//	else{
	//		pq_WP_overlap_array[i] = false;
	//	}
	//}

	size_t counter = 0;
	for (size_t idx=0; idx<Nq_WP*Np_WP*Nq_WP*Np_WP;idx++){
		if (pq_WP_overlap_array[idx]==false){
			counter += 1;
		}
	}
	double P123_sparsity = (double) counter/(Nq_WP*Nq_WP*Np_WP*Np_WP);
	printf("   - %.3f%% of P123-matrix violates momentum-conservation \n", 100*P123_sparsity );

	double*  sin_phi_array = new double [Nq_WP*Np_WP*Nphi];
	double*  cos_phi_array = new double [Nq_WP*Np_WP*Nphi];
	for (size_t idx=0; idx<Nq_WP*Np_WP*Nphi; idx++){
		sin_phi_array[idx] = sin(phi_array[idx]);
		cos_phi_array[idx] = cos(phi_array[idx]);
	}

	/* Preallocate array if we use dense format. Otherwise (i.e. sparse) start with
	 * some reasonable guess (usually less than a percent), and expand if required. */
	size_t P123_dense_dim    = Np_WP * Nq_WP * Nalpha;
	size_t P123_dense_dim_sq = P123_dense_dim * P123_dense_dim;

	/* Number of elements each thread can hold before writing to disk.
	 * Default is 1 GB memory per thread, as this is a (somewhat) safe minimum one typically
	 * finds on most computers today (the number of cores usually is equal to, or smaller than, the memory in GB) */
	size_t tread_buffer_size = std::pow(2,30)/(sizeof(double) + 2* sizeof(int));
	
	int P123_omp_num_threads = omp_get_max_threads();
	if (P123_omp_num_threads>Nq_WP){
		P123_omp_num_threads = Nq_WP;
	}

	/* Pointer-arrays for each OpenMP thread */
	double** P123_val_array_omp = NULL;
	int**    P123_row_array_omp = NULL;
	int**    P123_col_array_omp = NULL;
	size_t*  P123_dim_array_omp = NULL;
	/* Thread File Count (TFC) array. Used to keep track of which file the current
	 * thread will store to so as to empty contents (needed to prevent running out of memory) */
	int* 	 P123_TFC_array_omp	= NULL;

	if (use_dense_format){
		*P123_val_dense_array = new double [P123_dense_dim_sq];
		P123_dim = P123_dense_dim;
	}
	else{
		/* The sparse dimension is determined by counting, see the loops below */
		P123_dim = 0;

		P123_val_array_omp = new double* [P123_omp_num_threads];
		P123_row_array_omp = new int*    [P123_omp_num_threads];
		P123_col_array_omp = new int*    [P123_omp_num_threads];
		P123_dim_array_omp = new size_t  [P123_omp_num_threads];
		P123_TFC_array_omp = new int     [P123_omp_num_threads];

		for (int thread_idx=0; thread_idx<P123_omp_num_threads; thread_idx++){
			P123_val_array_omp[thread_idx] = new double [tread_buffer_size];
			P123_row_array_omp[thread_idx] = new int    [tread_buffer_size];
			P123_col_array_omp[thread_idx] = new int    [tread_buffer_size];
			P123_dim_array_omp[thread_idx] = P123_dim;
			P123_TFC_array_omp[thread_idx] = 0;
		}
	}

	//int Gtilde_subarray_size = Np_per_WP * Nq_per_WP * Nx_Gtilde;
	int Gtilde_subarray_size = Nphi * Nx;


	int row_step_length = P123_dense_dim/100;
	int row_multiplier = 1;

	/* Code segment to check sparse matrix size */
	double P123_density   = 1. - P123_sparsity;
	size_t mem_check_nnz  = P123_density * P123_dense_dim_sq;
	size_t mem_check_num_doubles = 2 * mem_check_nnz;
	size_t mem_check_num_ints    = 4 * mem_check_nnz; 
	double mem_check_size_in_GB  = (mem_check_num_doubles*sizeof(double) + mem_check_num_ints*sizeof(int))/std::pow(2.0,30);
	printf("   - Maximum memory required for P123 calculation + storage: %.2f GB \n", mem_check_size_in_GB);
	printf("     - Checking if required memory is available ... \n");
	double* mem_check_array_doubles = NULL;
	int*    mem_check_array_ints    = NULL;
	try{
		mem_check_array_doubles = new double [mem_check_num_doubles];
		mem_check_array_ints    = new int    [mem_check_num_ints];
	}
	catch (...) {
		raise_error("     - Memory allocation failed.");
	}
	delete [] mem_check_array_doubles;
	delete [] mem_check_array_ints;
	printf("     - Confirmed. \n");

	/* Start of P123 parallel calculation */
	printf("   - Initiating P123-matrix calculation ... \n");
	printf("     - Running OpenMP on %d threads \n", P123_omp_num_threads);
	fflush(stdout);
	auto timestamp_P123_tot_start = chrono::system_clock::now();
	
	int    num_rows_calculated  [P123_omp_num_threads];
	double tot_calculation_time [P123_omp_num_threads];
	for (int thread_idx=0; thread_idx<P123_omp_num_threads; thread_idx++){
		num_rows_calculated [thread_idx] = 0;
		tot_calculation_time[thread_idx] = 0;
	}
	
	#pragma omp parallel
	{
	int thread_idx = omp_get_thread_num();

	//#pragma omp for
	//for (int row=0; row<P123_dense_dim; row+=stplngth){
	//	for (int col=0; col<P123_dense_dim; col+=stplngth){
	//		double prob_nnz = (double) rand() / RAND_MAX;
	//		int P123_dim_omp = P123_dim_array_omp[thread_idx];;
	//		double P123_element = (double) rand() / RAND_MAX;
	//		(P123_val_array_omp[thread_idx])[P123_dim_omp] = P123_element;
	//		(P123_row_array_omp[thread_idx])[P123_dim_omp] = row;
	//		(P123_col_array_omp[thread_idx])[P123_dim_omp] = col;
	//		P123_dim_omp += 1;
	//		P123_dim_array_omp[thread_idx] = P123_dim_omp;
	//		if ( P123_dim_omp>=current_array_dim ){  // This should occur a small amount of the time
	//			increase_sparse_array_size(&P123_val_array_omp[thread_idx], current_array_dim, sparse_step_length);
	//			increase_sparse_array_size(&P123_row_array_omp[thread_idx], current_array_dim, sparse_step_length);
	//			increase_sparse_array_size(&P123_col_array_omp[thread_idx], current_array_dim, sparse_step_length);
	//			current_array_dim += sparse_step_length;
	//		}
	//		if (P123_dim_omp>50806056/16){
	//			break;
	//		}
	//	}
	//}
	//}

	int  	P123_row_idx = 0;
	int  	P123_col_idx = 0;
	size_t  pq_WP_idx	 = 0;
	bool    WP_overlap	 = false;
	double  P123_element = 0;

	int L_2N	   = 0;
	int L_1N	   = 0;
	int L_2N_prime = 0;
	int L_1N_prime = 0;
	
	double Gtilde_subarray [Gtilde_subarray_size];

	#pragma omp for
	/* <X_i'j'^alpha'| - loops (rows of P123) */
	for (int qp_idx_WP = 0; qp_idx_WP < Nq_WP; qp_idx_WP++){
		for (int pp_idx_WP = 0; pp_idx_WP < Np_WP; pp_idx_WP++){

			/* Progress printout by thread 0 */
			if (thread_idx==0){
				int num_rows_count = 0;
				for (int i=0; i<P123_omp_num_threads; i++){
					num_rows_count += num_rows_calculated[i];
				}

				auto timestamp_P123_tot_inter = chrono::system_clock::now();
				chrono::duration<double> time_P123_inter = timestamp_P123_tot_inter - timestamp_P123_tot_start;
				double P123_current_time = time_P123_inter.count();
				double P123_av_row_time  = P123_current_time/num_rows_count;
				double P123_est_completion_time = (P123_dense_dim-num_rows_count)*P123_av_row_time/3600.;
				
				printf("\r     - Calculated %d of %d rows. Av. time per row: %.3f s. Est. completion time: %.1f h", num_rows_count, P123_dense_dim, P123_av_row_time, P123_est_completion_time);
				fflush(stdout);
			}

			for (int alphap_idx = 0; alphap_idx < 1; alphap_idx++){
	
				/* |X_ij^alpha> - loops (columns of P123) */
				for (int alpha_idx = 0; alpha_idx < Nalpha; alpha_idx++){
	
					L_2N = L_2N_array[alphap_idx];
					L_1N = L_1N_array[alphap_idx];
	
					L_2N_prime = L_2N_array[alpha_idx];
					L_1N_prime = L_1N_array[alpha_idx];
					
					if (production_run){
						calculate_Gtilde_subarray_polar(Gtilde_subarray,
												  	    &Atilde_store[alphap_idx*Nalpha*(Lmax+1) + alpha_idx*(Lmax+1)],
								   					    Nx, x_array,
								   					    Nphi,
								   					    &sin_phi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi],
								   					    &cos_phi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi],
								   					    L_2N, L_2N_prime,
												  	    L_1N, L_1N_prime,
												  	    two_J_3N);
					}
	
					for (int q_idx_WP = 0; q_idx_WP < Nq_WP; q_idx_WP++){
						for (int p_idx_WP = 0; p_idx_WP < Np_WP; p_idx_WP++){
	
							/* Unique index for current combination of WPs */
							pq_WP_idx = qp_idx_WP*Np_WP*Nq_WP*Np_WP
									  + pp_idx_WP*Nq_WP*Np_WP
									  +  q_idx_WP*Np_WP
									  +  p_idx_WP;
							WP_overlap = pq_WP_overlap_array[pq_WP_idx];

							/* Only calculate P123 if there is WP bin-overlap in Heaviside functions */
							if (WP_overlap){
								if (production_run){
									P123_element = calculate_P123_element_in_WP_basis_mod (Gtilde_subarray,
											   												p_idx_WP, q_idx_WP,
											   												pp_idx_WP, qp_idx_WP,
											   												Np_WP,p_array_WP_bounds,
											   												Nq_WP,q_array_WP_bounds,
											   												Nx, x_array, wx_array,
											   												Nphi,
											   												&sin_phi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi],
											   												&cos_phi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi],
											   												&wphi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi]);
								}
								else{
									P123_row_idx = alphap_idx*Nq_WP*Np_WP + qp_idx_WP*Np_WP +  pp_idx_WP;
									P123_col_idx = alpha_idx*Nq_WP*Np_WP + q_idx_WP*Np_WP + p_idx_WP;

									P123_element = cos(P123_row_idx)*cos(P123_col_idx);
								}
								
								if (P123_element!=0){
									P123_row_idx = alphap_idx*Nq_WP*Np_WP + qp_idx_WP*Np_WP +  pp_idx_WP;
									P123_col_idx = alpha_idx*Nq_WP*Np_WP + q_idx_WP*Np_WP + p_idx_WP;

									if (use_dense_format){
										size_t P123_mat_idx = P123_row_idx*P123_dense_dim + P123_col_idx;
										(*P123_val_dense_array)[P123_mat_idx] = P123_element;
									}
									else{
										size_t P123_dim_omp = P123_dim_array_omp[thread_idx];

										(P123_val_array_omp[thread_idx])[P123_dim_omp] = P123_element;
										(P123_row_array_omp[thread_idx])[P123_dim_omp] = P123_row_idx;
										(P123_col_array_omp[thread_idx])[P123_dim_omp] = P123_col_idx;

										/* Increment sparse dimension (num of non-zero elements) */
										P123_dim_omp += 1;
										P123_dim_array_omp[thread_idx] = P123_dim_omp;

										/* Check if we have filled buffer array.
										 * If so, store array to disk and continue calculations.
										 * Each thread stores to their own file (filename includes thread index)
										 * such that there are no race hazard. Each write-to-disk creates a new file */
										/* Thread File Count (TFC) */
										int current_TFC = P123_TFC_array_omp[thread_idx];
										if ( P123_dim_omp>=tread_buffer_size ){

											std::string thread_filename = generate_subarray_file_name(two_J_3N, two_T_3N, P_3N,
																									  Np_WP, Nq_WP,
																									  J_2N_max,
																									  thread_idx,
																									  current_TFC);
											
											/* Store array */
											store_sparse_permutation_matrix_for_3N_channel_h5(P123_val_array_omp[thread_idx],
															  								  P123_row_array_omp[thread_idx],
															  								  P123_col_array_omp[thread_idx],
															  								  P123_dim_omp,
															  								  Np_WP, p_array_WP_bounds,
															  								  Nq_WP, q_array_WP_bounds,
															  								  Nalpha,
															  								  L_2N_array,
															  								  S_2N_array,
															  								  J_2N_array,
															  								  T_2N_array,
															  								  L_1N_array,
															  								  two_J_1N_array,
															  								  two_J_3N,
															  								  two_T_3N,
															  								  P_3N,
													   		  								  thread_filename,
																							  false);

											/* Re-set sparse-dimension. Old sparse-elements will be rewritten and/or not stored */
											P123_dim_array_omp[thread_idx] = 0;

											/* Increment the number of files stored for current thread */
											P123_TFC_array_omp[thread_idx] += 1;
										}
									}
								}
							}
						}
					}
				}

				num_rows_calculated[thread_idx] += 1;
			}
		}
	}
	}
	if (use_dense_format==false){
		/* Write remaining elements to file */
		for (int thread_idx=0; thread_idx<P123_omp_num_threads; thread_idx++){
			int current_TFC = P123_TFC_array_omp[thread_idx];

			std::string thread_filename = generate_subarray_file_name(two_J_3N, two_T_3N, P_3N,
																	  Np_WP, Nq_WP,
																	  J_2N_max,
																	  thread_idx,
																	  current_TFC);

			size_t num_elements_remaining = P123_dim_array_omp[thread_idx];

			/* Store array */
			store_sparse_permutation_matrix_for_3N_channel_h5(P123_val_array_omp[thread_idx],
							  								  P123_row_array_omp[thread_idx],
							  								  P123_col_array_omp[thread_idx],
							  								  num_elements_remaining,
							  								  Np_WP, p_array_WP_bounds,
							  								  Nq_WP, q_array_WP_bounds,
							  								  Nalpha,
							  								  L_2N_array,
							  								  S_2N_array,
							  								  J_2N_array,
							  								  T_2N_array,
							  								  L_1N_array,
							  								  two_J_1N_array,
							  								  two_J_3N,
							  								  two_T_3N,
							  								  P_3N,
					   		  								  thread_filename,
															  false);
			/* Increment the number of files stored for current thread */
			P123_TFC_array_omp[thread_idx] += 1;

			/* Deallocate thread arrays */
			delete [] P123_val_array_omp[thread_idx];
			delete [] P123_row_array_omp[thread_idx];
			delete [] P123_col_array_omp[thread_idx];
		}
		delete [] P123_val_array_omp;
		delete [] P123_row_array_omp;
		delete [] P123_col_array_omp;
		delete [] P123_dim_array_omp;
	}
	auto timestamp_P123_tot_end = chrono::system_clock::now();
	chrono::duration<double> time_P123_tot = timestamp_P123_tot_end - timestamp_P123_tot_start;
	double P123_tot_time = time_P123_tot.count();
	printf("\n     - Done. Tot. time used: %.1f h \n", P123_tot_time/3600.); fflush(stdout);

	/* Contract arrays to minimal size (number of non-zero elements) */
	printf("   - Merging P123 parallel-distributed arrays into a single COO-format structure ... \n"); fflush(stdout);
	if (use_dense_format==false){
		printf("   - Reading dimensions from parallel thread files \n"); fflush(stdout);
		/* Read dimensions from sparse subarray files */
		for (int thread_idx=0; thread_idx<P123_omp_num_threads; thread_idx++){
			int TFC_max = P123_TFC_array_omp[thread_idx];
			for (int current_TFC=0; current_TFC<TFC_max; current_TFC++){
				/* Generate filename for current thread_idx and TFC */
				std::string thread_filename = generate_subarray_file_name(two_J_3N, two_T_3N, P_3N,
																		  Np_WP, Nq_WP,
																		  J_2N_max,
																		  thread_idx,
																		  current_TFC);
				
				/* Convert std::string filename to char */
				char filename_char[300];
				std::strcpy(filename_char, thread_filename.c_str());
				
				/* Retrieve number of P123-elements in current file */
				unsigned long long current_P123_sparse_dim = 0;
				read_ULL_integer_from_h5(current_P123_sparse_dim, "P123_sparse_dim", filename_char);

				P123_dim += current_P123_sparse_dim;
			}
		}
		printf("   - Total number of non-zero elements: %zu \n", P123_dim); fflush(stdout);
		
		/* Allocate required memory to fite whole P123-matrix */
		double required_mem = P123_dim * (sizeof(double) + 2*sizeof(int))/std::pow(2.0,30);
		printf("   - Allocating necessary arrays (requires %.2f GB) \n", required_mem); fflush(stdout);
		try{
			*P123_val_sparse_array = new double [P123_dim];
			*P123_row_array        = new int    [P123_dim];
			*P123_col_array        = new int    [P123_dim];
		}
		catch (...) {
			raise_error("Failed. Memory exceeded in sparse memory allocation");
		}

		/* Read elements from sparse subarray files */
		printf("   - Reading elements (in random ordering) from parallel thread files \n"); fflush(stdout);
		size_t nnz_counter = 0;
		for (int thread_idx=0; thread_idx<P123_omp_num_threads; thread_idx++){
			int TFC_max = P123_TFC_array_omp[thread_idx];
			for (int current_TFC=0; current_TFC<TFC_max; current_TFC++){
				/* Generate filename for current thread_idx and TFC */
				std::string thread_filename = generate_subarray_file_name(two_J_3N, two_T_3N, P_3N,
																		  Np_WP, Nq_WP,
																		  J_2N_max,
																		  thread_idx,
																		  current_TFC);
				
				/* Retrieve P123-elements in current file and write to temporary arrays*/
				double* P123_sparse_val_subarray = NULL;
				int* 	P123_sparse_row_subarray = NULL;
				int* 	P123_sparse_col_subarray = NULL;
				size_t  current_P123_sparse_dim  = 0;
				read_sparse_permutation_matrix_for_3N_channel_h5(&P123_sparse_val_subarray,
															     &P123_sparse_row_subarray,
															     &P123_sparse_col_subarray,
															     current_P123_sparse_dim,
															     Np_WP, p_array_WP_bounds,
															     Nq_WP, q_array_WP_bounds,
															     Nalpha,
															     L_2N_array,
															     S_2N_array,
															     J_2N_array,
															     T_2N_array,
															     L_1N_array,
															     two_J_1N_array,
															     two_J_3N,
															     two_T_3N,
															     P_3N,
													   		     thread_filename,
																 false);

				/* Append RANDOM ORDER elements to input arrays */
				for (size_t idx=0; idx<current_P123_sparse_dim; idx++){
					(*P123_val_sparse_array)[idx+nnz_counter] = P123_sparse_val_subarray[idx];
					(*P123_row_array)       [idx+nnz_counter] = P123_sparse_row_subarray[idx];
					(*P123_col_array)       [idx+nnz_counter] = P123_sparse_col_subarray[idx];
				}

				nnz_counter += current_P123_sparse_dim;

				delete [] P123_sparse_val_subarray;
				delete [] P123_sparse_row_subarray;
				delete [] P123_sparse_col_subarray;
			}
		}

		/* Index array to sort */
		printf("   - Creating sorting vector for row-major COO-format \n"); fflush(stdout);
		std::vector<size_t> P123_idx_vector (P123_dim, 0);
		for (size_t idx=0; idx<P123_dim; idx++){
			P123_idx_vector[idx] = (size_t)(*P123_row_array)[idx]*P123_dense_dim + (size_t)(*P123_col_array)[idx];
		}

		/* Sort indices using template (Warning - lambda expression - requires C++11 or newer compiler) */
		auto sorted_indices = sort_indexes(P123_idx_vector);

		///* Sorting test */
		//std::vector<int> test_vec (10, 0);;
		//for (int idx=0; idx<10; idx++){
		//	test_vec[idx] = -(5-idx);
		//	std::cout << "test_vec: " << test_vec[idx] << std::endl;
		//}
		//auto sorted_test = sort_indexes(test_vec);
		//for (int idx=0; idx<10; idx++){
		//	std::cout << "test_vec: " << test_vec[sorted_test[idx]] << std::endl;
		//}

		/* Sort row indices */
		printf("     - Sorting random-order row indices \n"); fflush(stdout);
		int* sparse_idx_array_temp = new int [P123_dim];
		for (size_t i=0; i<P123_dim; i++){
			sparse_idx_array_temp[i] = (*P123_row_array)[sorted_indices[i]];
		}
		std::copy(sparse_idx_array_temp, sparse_idx_array_temp + P123_dim, *P123_row_array);

		//for (int i=0; i<P123_dim-1; i++){
		//	if (sparse_idx_array_temp[i]>sparse_idx_array_temp[i+1]){
		//		std::cout << i << " " << sparse_idx_array_temp[i] << " " << sparse_idx_array_temp[i+1] << std::endl;
		//		std::cout << i << " " << (*P123_row_array)[sorted_indices[i]] << " " << (*P123_row_array)[sorted_indices[i+1]] << std::endl;
		//		std::cout << i << " " << (*P123_row_array)[i] << " " << (*P123_row_array)[i+1] << std::endl;
		//		raise_error("row index error");
		//	}
		//}

		/* Sort column indices */
		printf("     - Sorting random-order column indices \n"); fflush(stdout);
		for (size_t i=0; i<P123_dim; i++){
			sparse_idx_array_temp[i] = (*P123_col_array)[sorted_indices[i]];
		}
		std::copy(sparse_idx_array_temp, sparse_idx_array_temp + P123_dim, *P123_col_array);

		delete [] sparse_idx_array_temp;
		double* sparse_val_array_temp = new double [P123_dim];

		/* Sort values */
		printf("     - Sorting random-order P123-values \n"); fflush(stdout);
		for (size_t i=0; i<P123_dim; i++){
			sparse_val_array_temp[i] = (*P123_val_sparse_array)[sorted_indices[i]];
		}
		std::copy(sparse_val_array_temp, sparse_val_array_temp + P123_dim, *P123_val_sparse_array);
		delete [] sparse_val_array_temp;
		printf("     - Done \n"); fflush(stdout);
	}

	/* Delete all temporary arrays */
	delete [] pq_WP_overlap_array;
	delete [] Atilde_store;
	delete [] SixJ_array;
}

void calculate_permutation_matrices_for_all_3N_channels(double** P123_sparse_val_array,
														int**    P123_sparse_row_array,
														int**    P123_sparse_col_array,
														size_t&  P123_sparse_dim,
														bool     production_run,
														int      Np_WP, double *p_array_WP_bounds,
														int      Nq_WP, double *q_array_WP_bounds,
														int      Nx, double* x_array, double* wx_array,
														int      Nphi,
														int      J_2N_max,
														int      Nalpha,
														int*     L_2N_array,
														int*     S_2N_array,
														int*     J_2N_array,
														int*     T_2N_array,
														int*     L_1N_array,
														int*     two_J_1N_array,
														int      two_J_3N,
														int      two_T_3N,
														int      P_3N){
	
	bool print_content = true;

	/* Either use a dense-block or allocate space dynamically
	 * as non-zero elements are found. The 2nd approach is a bit slower,
	 * but necessary due to the typical size of the dense-block being unhandleable. */
	bool use_dense_format = false;
	
	/* Dense and sparse array dimensions */
	size_t P123_array_dense_dim  = Nalpha * Np_WP * Nq_WP;
	size_t P123_array_dim        = 0;

	/* Temporary dense matrix to be filled in subroutine if use_dense_format=true */
	double* P123_val_dense_array = NULL;

	if (print_content){
		printf(" - Begin calculating 3N-channel permutation matrix \n");
	}
	calculate_permutation_matrix_for_3N_channel(&P123_val_dense_array,
												P123_sparse_val_array,
												P123_sparse_row_array,
												P123_sparse_col_array,
												P123_array_dim,
												use_dense_format,
												production_run,
							  					Np_WP, p_array_WP_bounds,
							  					Nq_WP, q_array_WP_bounds,
							  					Nx, x_array, wx_array,
												Nphi,
												J_2N_max,
							  					Nalpha,
							  					L_2N_array,
							  					S_2N_array,
							  					J_2N_array,
							  					T_2N_array,
							  					L_1N_array,
							  					two_J_1N_array,
												two_J_3N,
												two_T_3N,
												P_3N);

	if (print_content){
		printf("   - Calculation finished. \n");
		fflush(stdout);
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
													P123_sparse_dim);
	}
	else{
		P123_sparse_dim = P123_array_dim;
	}
		
	if (print_content){
		size_t P123_array_dense_dim_sq = P123_array_dense_dim*P123_array_dense_dim;
		double P_123_subarray_density = (double) P123_sparse_dim / P123_array_dense_dim_sq;

		printf(" - Dense dimension:   %dx%d \n", P123_array_dense_dim, P123_array_dense_dim);
		printf(" - Non-zero elements: %zu \n",   P123_sparse_dim);
		printf(" - Density:           %.3f %% \n",  100*P_123_subarray_density);
		fflush(stdout);
	}
		
	delete [] P123_val_dense_array;
}
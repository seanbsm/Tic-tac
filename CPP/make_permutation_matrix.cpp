
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

int factorial(int n){
    if(n > 1)
        return n * factorial(n - 1);
    else
        return 1;
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

//void calculate_Gtilde_subarray_cart(double* Gtilde_subarray,
//									double* Atilde_subarray,
//									int  Nx, double* x_array,
//									int  Nq, double* q_array,
//									int  Np, double* p_array,
//									int  L_2N, int L_2N_prime,
//									int  L_1N, int L_1N_prime,
//									int  two_J_3N){
//
//	bool print_Gtilde_progress = false;
//
//	long int fullsize = Np*Nq*Nx;
//	long int counter = 0;
//	int frac_n, frac_o=0;
//	
//	for (MKL_INT64 p_idx=0; p_idx<Np; p_idx++){
//		for (MKL_INT64 q_idx=0; q_idx<Nq; q_idx++){
//			for (MKL_INT64 x_idx=0; x_idx<Nx; x_idx++){
//				Gtilde_subarray[p_idx*Nq*Nx + q_idx*Nx + x_idx]
//					= Gtilde_subarray_new (p_array[p_idx],
//										   q_array[q_idx],
//										   x_array[x_idx],
//										   L_2N, L_2N_prime,
//										   L_1N, L_1N_prime,
//										   Atilde_subarray,
//										   two_J_3N);
//
//				counter += 1;
//			
//				if (print_Gtilde_progress){
//					frac_n = (100*counter)/fullsize;
//					if (frac_n>frac_o){std::cout << frac_n << "%" << std::endl; frac_o=frac_n;}
//				}
//			}
//		}
//	}
//}

void calculate_Gtilde_subarray_polar(double* Gtilde_subarray,
								     double* Atilde_subarray,
								     int  Nx, double* x_array,
								     int  Nphi,
								     double* sin_phi_subarray,
								     double* cos_phi_subarray,
								     int  L_2N, int L_2N_prime, int max_L_2N,
								     int  L_1N, int L_1N_prime, int max_L_1N,
								     int  two_J_3N,
									 double* ClebschGordan_data,
									 double* gsl_Plm_1_subarray, size_t gsl_Plm_1_stplen,
									 double* gsl_Plm_2_subarray, size_t gsl_Plm_2_stplen,
									 double* gsl_Plm_3_subarray, size_t gsl_Plm_3_stplen,
									 double* prefac_L_array,
									 double* prefac_l_array,
									 int two_jmax_Clebsch){

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
									   L_2N, L_2N_prime, max_L_2N,
									   L_1N, L_1N_prime, max_L_1N,
									   Atilde_subarray,
									   two_J_3N,
									   ClebschGordan_data,
									   &gsl_Plm_1_subarray[x_idx*gsl_Plm_1_stplen],
									   &gsl_Plm_2_subarray[(phi_idx*Nx + x_idx)*gsl_Plm_2_stplen],
									   &gsl_Plm_3_subarray[(phi_idx*Nx + x_idx)*gsl_Plm_3_stplen],
									   prefac_L_array,
									   prefac_l_array,
									   two_jmax_Clebsch);

			counter += 1;
			
			if (print_Gtilde_progress){
				frac_n = (100*counter)/fullsize;
				if (frac_n>frac_o){std::cout << frac_n << "%" << std::endl; frac_o=frac_n;}
			}
		}
	}
}

std::string generate_subarray_file_name(int two_J_3N, int P_3N,
										int Np_WP, int Nq_WP,
										int J_2N_max,
										int thread_idx,
										int current_TFC,
										std::string P123_folder){

	std::string filename = P123_folder + "/P123_subsparse_JP_"
						 + to_string(two_J_3N) + "_" + to_string(P_3N)
						 + "_Np_" + to_string(Np_WP) + "_Nq_" + to_string(Nq_WP)
						 + "_J2max_" + to_string(J_2N_max) + "_TFC_" + to_string(thread_idx)
						 + "_" + to_string(current_TFC) + ".h5";

	return filename;
}

void calculate_permutation_elements_for_3N_channel(double** P123_val_dense_array,
												   int*		max_TFC_array,
												   bool     use_dense_format,
												   bool     production_run,
												   int      Np_WP, double *p_array_WP_bounds,
												   int      Nq_WP, double *q_array_WP_bounds,
												   int      Nx, double* x_array, double* wx_array,
												   int      Nphi,
												   int      J_2N_max,
												   pw_3N_statespace pw_states,
												   run_params run_parameters,
												   std::string P123_folder){
	
	/* Make local pointers & variables */
	int  Nalpha			= pw_states.Nalpha;
	int* L_2N_array		= pw_states.L_2N_array;
	int* S_2N_array		= pw_states.S_2N_array;
	int* J_2N_array		= pw_states.J_2N_array;
	int* T_2N_array		= pw_states.T_2N_array;
	int* L_1N_array		= pw_states.L_1N_array;
	int* two_J_1N_array = pw_states.two_J_1N_array;
	int* two_T_3N_array	= pw_states.two_T_3N_array;
	int  two_J_3N    	= pw_states.two_J_3N_array[0];
	int  P_3N   	    = pw_states.P_3N_array[0];
	
	bool print_content = true;

	/* Notation change */
	int two_J = two_J_3N;
	//int two_T = two_T_3N;
	
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
	
	//int l_interpolate_max = l_interpolate_max = 2 * (lmax - 3) + 3;
	//if (print_content){
	//	std::cout << "   - lmax = " << lmax << ", l_interpolate_max = " << l_interpolate_max << "\n";
	//}

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
	if (print_content){std::cout << "   - lmax      = " <<  lmax << std::endl;}

	if (print_content){
		printf("   - Precalculating Wigner 6j-symbols \n");
	}
	int SixJ_size = (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1);
	if (SixJ_size < 0){
		raise_error("SixJ_array in make_permutation_matrix had negative size, likely an integer overflow. Check your dimensions.");
	}
	printf("     - Total prestore requirement is %zu doubles. Allocating arrays ... \n", SixJ_size);
	double *SixJ_array = new double [SixJ_size];
	printf("     - Success. Calculating ... \n");
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
	printf("     - Done \n");

	if (print_content){
		printf("   - Precalculating Atilde \n");
	}
	MKL_INT64 Atilde_N = Jj_dim * Jj_dim * (Lmax + 1);
	printf("     - Total prestore requirement is %zu doubles. Allocating arrays ... \n", Atilde_N);
	double *Atilde_store = new double[Atilde_N];
	printf("     - Success. Calculating ... \n");
	for (MKL_INT64 i=0; i<Atilde_N; i++){
		Atilde_store[i] = 0.0;
	}
	for (int alpha = 0; alpha <= Jj_dim - 1; alpha++){
		int two_T_3N_alpha = two_T_3N_array[alpha];
		for (int alphaprime = 0; alphaprime <= Jj_dim - 1; alphaprime++){
			int two_T_3N_alphaprime = two_T_3N_array[alphaprime];
			if (two_T_3N_alpha==two_T_3N_alphaprime){
				int two_T = two_T_3N_alpha;
				for (int Ltotal = 0; Ltotal <= Lmax; Ltotal++){
					Atilde_store[alpha * Jj_dim * (Lmax + 1) + alphaprime * (Lmax + 1) + Ltotal] = Atilde (alpha, alphaprime, Ltotal, Jj_dim, L12_Jj, l3_Jj, J12_Jj, two_j3_Jj, S12_Jj, T12_Jj, two_J, two_T, SixJ_array, two_jmax_SixJ);
				}
			}
		}
	}
	if (print_content){
		printf("     - Done \n");
	}

	/* Precalculate Clebsch-Gordan coefficients in geometric function */
	int two_jmax_Clebsch = 2 * lmax;
    int jmax_Clebsch = lmax;

    // prestore Clebsch Gordan coefficients
    printf("   - Precalculating Clebsch-Gordan coefficients \n");
	size_t  ClebschGordan_size = (two_jmax_Clebsch + 1) * (two_jmax_Clebsch + 1) * (two_jmax_Clebsch + 1) * (two_jmax_Clebsch + 1) * (2 * two_jmax_Clebsch + 1);
    printf("     - Total prestore requirement is %zu doubles. Allocating arrays ... \n", ClebschGordan_size);
	double* ClebschGordan_data = new double [ClebschGordan_size];
	printf("     - Success. Calculating ... \n");
    #pragma omp parallel
    {
        #pragma omp for collapse(3)

        for (int two_j1 = 0; two_j1 <= two_jmax_Clebsch; two_j1++)
        {
            for (int two_j2 = 0; two_j2 <= two_jmax_Clebsch; two_j2++)
            {
                for (int two_j3 = 0; two_j3 <= 2 * two_jmax_Clebsch; two_j3++) // tabulate these l up to 2*lmax
                {
                    for (int two_m1 = -two_j1; two_m1 <= two_j1; two_m1++)
                    {
                        for (int two_m2 = -two_j2; two_m2 <= two_j2; two_m2++)
                        {
                            ClebschGordan_data[((two_j1 + 1) * (two_j1 + 1) + two_m1 - two_j1 - 1) * (two_jmax_Clebsch + 1) * (two_jmax_Clebsch + 1) * (2 * two_jmax_Clebsch + 1)
                                               + ((two_j2 + 1) * (two_j2 + 1) + two_m2 - two_j2 - 1) * (2 * two_jmax_Clebsch + 1)
                                               + two_j3]
                                = ClebschGordan(two_j1, two_j2, two_j3, two_m1, two_m2, two_m1 + two_m2);
                        }
                    }
                }
            }
        }
    }
	printf("     - Done \n");

	/* END OF OLD CODE SEGMENT WITH OLD VARIABLE-NOTATION */

	/* Precalculate overlapping bins and where p_bar and q_bar are non-zero */
	size_t	 num_WP_cells		 = (size_t)Nq_WP*Nq_WP*Np_WP*Np_WP;
	bool*    pq_WP_overlap_array = new bool   [num_WP_cells];
	double*   phi_array      	 = new double [(size_t)Nq_WP*Np_WP*Nphi];
	double*  wphi_array      	 = new double [(size_t)Nq_WP*Np_WP*Nphi];
	
	printf("   - Precalculating momentum conservations \n");
	fflush(stdout);
	if (true){
	#pragma omp parallel
	{
		#pragma omp for
	for (size_t qp_idx_WP=0; qp_idx_WP<Nq_WP; qp_idx_WP++){
		double qp_l = q_array_WP_bounds[qp_idx_WP];
		double qp_u = q_array_WP_bounds[qp_idx_WP+1];
		for (size_t pp_idx_WP=0; pp_idx_WP<Np_WP; pp_idx_WP++){
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
			for (size_t q_idx_WP=0; q_idx_WP<Nq_WP; q_idx_WP++){
				double q_l = q_array_WP_bounds[q_idx_WP];
				double q_u = q_array_WP_bounds[q_idx_WP+1];
				for (size_t p_idx_WP=0; p_idx_WP<Np_WP; p_idx_WP++){
					double p_l = p_array_WP_bounds[p_idx_WP];
					double p_u = p_array_WP_bounds[p_idx_WP+1];
	
					bool WP_overlap = false;
					/* Ensure possible phi boundaries */
					int hit_counter = 0;
					if (phi_lower<phi_upper){
						/* Search for on-shell elements */
						for (size_t phi_idx=0; phi_idx<Nphi; phi_idx++){
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
					size_t step_length_1 = (size_t) Np_WP*Nq_WP*Np_WP;
					size_t step_length_2 = (size_t) 	  Nq_WP*Np_WP;
					size_t step_length_3 = (size_t)		        Np_WP;
					size_t pq_WP_idx = (size_t) qp_idx_WP*step_length_1
							  	  	 		  + pp_idx_WP*step_length_2
							  	  	 		  +  q_idx_WP*step_length_3
							  	  	 		  +  p_idx_WP;
					//std::cout << pq_WP_idx << " " << qp_idx_WP << " " << pp_idx_WP << " " << q_idx_WP << " " << p_idx_WP << std::endl;
					pq_WP_overlap_array[pq_WP_idx] = WP_overlap;
				}
			}
		}
	}
	}
	}
	else{
		#pragma omp parallel
		{
		#pragma omp for
		for (size_t qp_idx_WP=0; qp_idx_WP<Nq_WP; qp_idx_WP++){
			double qp_l = q_array_WP_bounds[qp_idx_WP];
			double qp_u = q_array_WP_bounds[qp_idx_WP+1];
			for (size_t pp_idx_WP=0; pp_idx_WP<Np_WP; pp_idx_WP++){
				double pp_l = p_array_WP_bounds[pp_idx_WP];
				double pp_u = p_array_WP_bounds[pp_idx_WP+1];

				double phi_lower = atan(pp_l/qp_u);
				double phi_upper = atan(pp_u/qp_l);

				/* Create phi-mesh */
				calc_gauss_points ( &phi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi],
								   &wphi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi],
								   phi_lower, phi_upper,
								   Nphi);

				for (size_t q_idx_WP=0; q_idx_WP<Nq_WP; q_idx_WP++){
					double q_l = q_array_WP_bounds[q_idx_WP];
					double q_u = q_array_WP_bounds[q_idx_WP+1];
					for (size_t p_idx_WP=0; p_idx_WP<Np_WP; p_idx_WP++){
						double p_l = p_array_WP_bounds[p_idx_WP];
						double p_u = p_array_WP_bounds[p_idx_WP+1];
	
						double pi1_min = pi1_tilde(pp_l, qp_l, -1.0);
						double pi1_max = pi1_tilde(pp_u, qp_u, +1.0);
						bool p_in_pi1 = ( (pi1_min<=p_l && p_l<=pi1_max) || (pi1_min<=p_u && p_u<=pi1_max) );

						double pi2_min = pi2_tilde(pp_l, qp_l, +1.0);
						double pi2_max = pi2_tilde(pp_u, qp_u, -1.0);
						bool q_in_pi2 = ( (pi2_min<=q_l && q_l<=pi2_max) || (pi2_min<=q_u && q_u<=pi2_max) );

						bool WP_overlap = false;
						if ( p_in_pi1 && q_in_pi2 ){
							WP_overlap = true;
						}
					
						/* Unique index for current combination of WPs */
						size_t step_length_1 = (size_t) Np_WP*Nq_WP*Np_WP;
						size_t step_length_2 = (size_t) 	  Nq_WP*Np_WP;
						size_t step_length_3 = (size_t)		        Np_WP;
						size_t pq_WP_idx = (size_t) qp_idx_WP*step_length_1
								  	  	 		  + pp_idx_WP*step_length_2
								  	  	 		  +  q_idx_WP*step_length_3
								  	  	 		  +  p_idx_WP;
						//std::cout << pq_WP_idx << " " << qp_idx_WP << " " << pp_idx_WP << " " << q_idx_WP << " " << p_idx_WP << std::endl;
						pq_WP_overlap_array[pq_WP_idx] = WP_overlap;
					}
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
	for (size_t idx=0; idx<num_WP_cells;idx++){
		if (pq_WP_overlap_array[idx]==false){
			counter += 1;
		}
	}
	double P123_sparsity = (double) counter/num_WP_cells;
	printf("     - %.3f%% of P123-matrix violates momentum-conservation \n", 100*P123_sparsity );

	double*  sin_phi_array = new double [Nq_WP*Np_WP*Nphi];
	double*  cos_phi_array = new double [Nq_WP*Np_WP*Nphi];
	for (size_t idx=0; idx<Nq_WP*Np_WP*Nphi; idx++){
		sin_phi_array[idx] = sin(phi_array[idx]);
		cos_phi_array[idx] = cos(phi_array[idx]);
	}


	/* Precalculate Legendre polynomials in geometric function */
	printf("   - Precalculating Legendre polynomials \n");
	fflush(stdout);
	size_t len_Plm_array_L = gsl_sf_legendre_array_n(max_L12);
	size_t len_Plm_array_l = gsl_sf_legendre_array_n(max_l3);
	printf("     - L12-array requires %zu doubles per (phi,x). \n", len_Plm_array_L);
	printf("     - l3-array  requires %zu doubles per (phi,x). \n", len_Plm_array_l);

	size_t tot_size_needed = (len_Plm_array_L+len_Plm_array_l)*Np_WP*Nq_WP*Nphi*Nx+len_Plm_array_l*Nx;
	printf("     - Total prestore requirement is %zu doubles. Allocating arrays \n", tot_size_needed);
	fflush(stdout);

	double* gsl_Plm_1_array = new double [len_Plm_array_l * Nx];
	double* gsl_Plm_2_array = new double [len_Plm_array_L * Np_WP * Nq_WP * Nphi * Nx];
	double* gsl_Plm_3_array = new double [len_Plm_array_l * Np_WP * Nq_WP * Nphi * Nx];

	size_t gsl_Plm_1_stplen = len_Plm_array_l;
	size_t gsl_Plm_2_stplen = len_Plm_array_L;
	size_t gsl_Plm_3_stplen = len_Plm_array_l;

	/* Calculate all Plm needed for given input angles using gsl.
	 * The 1st argument is set to GSL_SF_LEGENDRE_NONE, meaning we calculate the "unnormalized associated Legendre polynomials Plm"
	 * The 4th argument is set to -1, which means we use the Condon-Shortley phase factor (-1)^m */
	printf("     - Calculating polynomials \n", tot_size_needed);
	fflush(stdout);
	for (int x_idx=0; x_idx<Nx; x_idx++){
		double x = x_array[x_idx];
		size_t idx_subarray_1 = x_idx;
		gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM, max_l3, x, -1, &gsl_Plm_1_array[idx_subarray_1 * gsl_Plm_1_stplen]);
	}
	for (size_t qp_idx_WP=0; qp_idx_WP<Nq_WP; qp_idx_WP++){
		for (size_t pp_idx_WP=0; pp_idx_WP<Np_WP; pp_idx_WP++){
			for (size_t phi_idx=0; phi_idx<Nphi; phi_idx++){
				double phi = phi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi + phi_idx];
				double sin_phi = std::sin(phi);
				double cos_phi = std::cos(phi);
				for (int x_idx=0; x_idx<Nx; x_idx++){
					double x = x_array[x_idx];
					double p = sin_phi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi + phi_idx];
					double q = cos_phi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi + phi_idx];
					double pi1 = pi1_tilde(p, q, x);
					double pi2 = pi2_tilde(p, q, x);
					double costheta1 = -(0.5 * p + 0.75 * q * x) / pi1;
					double costheta2 = (p - 0.5 * q * x) / pi2;

					/* Prevent numerical error in Plm */
					if ( costheta1>1 ){
						costheta1 = 1;
					}
					else if( costheta1<-1 ){
						costheta1 = -1;
					}
	
					/* Prevent numerical error in Plm */
					if ( costheta2>1 ){
						costheta2 = 1;
					}
					else if( costheta2<-1 ){
						costheta2 = -1;
					}

					size_t idx_subarray_2 = qp_idx_WP*Np_WP*Nphi*Nx + pp_idx_WP*Nphi*Nx + phi_idx*Nx + x_idx;
					size_t idx_subarray_3 = qp_idx_WP*Np_WP*Nphi*Nx + pp_idx_WP*Nphi*Nx + phi_idx*Nx + x_idx;
					gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM, max_L12, costheta1, -1, &gsl_Plm_2_array[idx_subarray_2 * gsl_Plm_2_stplen]);
					gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM, max_l3,  costheta2, -1, &gsl_Plm_3_array[idx_subarray_3 * gsl_Plm_3_stplen]);
				}
			}
		}
	}
	printf("     - Done \n");
	fflush(stdout);

	/* Arrays to hold prefactors to calculate Plm when m<0 (not calculated automatically by gsl)
	 * Use same indexing as gsl Plm-arrays */
	double* prefac_L_array  = new double [len_Plm_array_L];
	double* prefac_l_array  = new double [len_Plm_array_l];
	for (int L=0; L<max_L12+1; L++){
		for (int m=0; m<L+1; m++){
			int gsl_idx = gsl_sf_legendre_array_index(L, m);
			prefac_L_array[gsl_idx] = gsl_sf_pow_int(-1.0, m);//*(double)factorial(L-m)/factorial(L+m);
		}
	}
	for (int l=0; l<max_l3+1; l++){
		for (int m=0; m<l+1; m++){
			int gsl_idx = gsl_sf_legendre_array_index(l, m);
			prefac_l_array[gsl_idx] = gsl_sf_pow_int(-1.0, m);//*(double)factorial(l-m)/factorial(l+m);
		}
	}


	/* Preallocate array if we use dense format. Otherwise (i.e. sparse) start with
	 * some reasonable guess (usually less than a percent), and expand if required. */
	size_t P123_dense_dim    = Np_WP * Nq_WP * Nalpha;
	size_t P123_dense_dim_sq = P123_dense_dim * P123_dense_dim;

	/* Number of elements each thread can hold before writing to disk.
	 * Default is 1 GB memory per thread, as this is a (somewhat) safe minimum one typically
	 * finds on most computers today (the number of cores usually is equal to, or smaller than, the memory in GB) */
	size_t tread_buffer_size = std::pow(2,30)/(sizeof(double) + 2* sizeof(int));
	
	int P123_omp_num_threads = run_parameters.P123_omp_num_threads;

	/* Pointer-arrays for each OpenMP thread */
	double** P123_val_array_omp = NULL;
	int**    P123_row_array_omp = NULL;
	int**    P123_col_array_omp = NULL;
	size_t*  P123_dim_array_omp = NULL;

	if (use_dense_format){
		*P123_val_dense_array = new double [P123_dense_dim_sq];
	}
	else{
		P123_val_array_omp = new double* [P123_omp_num_threads];
		P123_row_array_omp = new int*    [P123_omp_num_threads];
		P123_col_array_omp = new int*    [P123_omp_num_threads];
		P123_dim_array_omp = new size_t  [P123_omp_num_threads];

		for (int thread_idx=0; thread_idx<P123_omp_num_threads; thread_idx++){
			P123_val_array_omp[thread_idx] = new double [tread_buffer_size];
			P123_row_array_omp[thread_idx] = new int    [tread_buffer_size];
			P123_col_array_omp[thread_idx] = new int    [tread_buffer_size];
			P123_dim_array_omp[thread_idx] = 0;

			max_TFC_array[thread_idx] = 0;
		}
	}

	//int Gtilde_subarray_size = Np_per_WP * Nq_per_WP * Nx_Gtilde;
	size_t Gtilde_subarray_size = Nphi * Nx;


	//int row_step_length = P123_dense_dim/100;
	//int row_multiplier = 1;

	/* Code segment to check sparse matrix size */
	double P123_density   = 1. - P123_sparsity;
	size_t mem_check_nnz  = P123_density * P123_dense_dim_sq;
	size_t mem_check_num_doubles = 2 * mem_check_nnz;
	size_t mem_check_num_ints    = 4 * mem_check_nnz; 
	double mem_check_size_in_GB  = (mem_check_num_doubles*sizeof(double) + mem_check_num_ints*sizeof(int))/std::pow(2.0,30);
	printf("   - Maximum memory required for P123 calculation + storage routines: %.2f GB \n", mem_check_size_in_GB);
	printf("     - Checking if required memory is available ... \n");
	double* mem_check_array_doubles = NULL;
	int*    mem_check_array_ints    = NULL;
	//try{
	//	mem_check_array_doubles = new double [mem_check_num_doubles];
	//	mem_check_array_ints    = new int    [mem_check_num_ints];
	//}
	//catch (...) {
	//	raise_error("     - Memory allocation failed.");
	//}
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
	
	omp_set_num_threads(P123_omp_num_threads);
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

	size_t	P123_row_idx = 0;
	size_t	P123_col_idx = 0;
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
	for (size_t qp_idx_WP = 0; qp_idx_WP < Nq_WP; qp_idx_WP++){
		for (size_t pp_idx_WP = 0; pp_idx_WP < Np_WP; pp_idx_WP++){

			/* Progress printout by thread 0 */
			if (thread_idx==0){
				int num_rows_count = 0;
				for (size_t i=0; i<P123_omp_num_threads; i++){
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

			for (size_t alphap_idx = 0; alphap_idx < Nalpha; alphap_idx++){
	
				/* |X_ij^alpha> - loops (columns of P123) */
				for (size_t alpha_idx = 0; alpha_idx < Nalpha; alpha_idx++){
	
					L_2N = L_2N_array[alphap_idx];
					L_1N = L_1N_array[alphap_idx];
	
					L_2N_prime = L_2N_array[alpha_idx];
					L_1N_prime = L_1N_array[alpha_idx];

					// ONLY USE S-WAVE (handy for Malfliet-Tjon debugging)
					//if (L_2N!=0 || L_2N_prime!=0){
					//	continue;
					//}
					
					if (production_run){
						calculate_Gtilde_subarray_polar(Gtilde_subarray,
												  	    &Atilde_store[alphap_idx*Nalpha*(Lmax+1) + alpha_idx*(Lmax+1)],
								   					    Nx, x_array,
								   					    Nphi,
								   					    &sin_phi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi],
								   					    &cos_phi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi],
								   					    L_2N, L_2N_prime, max_L12,
												  	    L_1N, L_1N_prime, max_l3,
												  	    two_J_3N,
														ClebschGordan_data,
														 gsl_Plm_1_array, gsl_Plm_1_stplen,
														&gsl_Plm_2_array[(qp_idx_WP*Np_WP*Nphi*Nx + pp_idx_WP*Nphi*Nx)*gsl_Plm_2_stplen], gsl_Plm_2_stplen,
														&gsl_Plm_3_array[(qp_idx_WP*Np_WP*Nphi*Nx + pp_idx_WP*Nphi*Nx)*gsl_Plm_3_stplen], gsl_Plm_3_stplen,
														prefac_L_array,
														prefac_l_array,
														two_jmax_Clebsch);
					}
	
					for (size_t q_idx_WP = 0; q_idx_WP < Nq_WP; q_idx_WP++){
						for (size_t p_idx_WP = 0; p_idx_WP < Np_WP; p_idx_WP++){
	
							/* Unique index for current combination of WPs */
							size_t step_length_1 = (size_t) Np_WP*Nq_WP*Np_WP;
							size_t step_length_2 = (size_t) 	  Nq_WP*Np_WP;
							size_t step_length_3 = (size_t)		        Np_WP;
							size_t pq_WP_idx = (size_t) qp_idx_WP*step_length_1
									  	  	 		  + pp_idx_WP*step_length_2
									  	  	 		  +  q_idx_WP*step_length_3
									  	  	 		  +  p_idx_WP;
							WP_overlap = pq_WP_overlap_array[pq_WP_idx];

							/* Only calculate P123 if there is WP bin-overlap in Heaviside functions */
							if (WP_overlap){
								if (production_run){
									P123_element = calculate_P123_element_in_WP_basis_mod (Gtilde_subarray,
											   												(int) p_idx_WP,  (int) q_idx_WP,
											   												(int) pp_idx_WP, (int) qp_idx_WP,
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
										int current_TFC = max_TFC_array[thread_idx];
										if ( P123_dim_omp>=tread_buffer_size ){

											std::string thread_filename = generate_subarray_file_name(two_J_3N, P_3N,
																									  Np_WP, Nq_WP,
																									  J_2N_max,
																									  thread_idx,
																									  current_TFC,
																									  P123_folder);
											
											/* Store array */
											store_sparse_permutation_matrix_for_3N_channel_h5(P123_val_array_omp[thread_idx],
															  								  P123_row_array_omp[thread_idx],
															  								  P123_col_array_omp[thread_idx],
															  								  P123_dim_omp,
															  								  Np_WP, p_array_WP_bounds,
															  								  Nq_WP, q_array_WP_bounds,
															  								  pw_states,
													   		  								  thread_filename,
																							  false);

											/* Re-set sparse-dimension. Old sparse-elements will be rewritten and/or not stored */
											P123_dim_array_omp[thread_idx] = 0;

											/* Increment the number of files stored for current thread */
											max_TFC_array[thread_idx] += 1;
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
	/* This changes back the maximal number of threads permitted by omp */
	omp_set_num_threads(omp_get_max_threads());

	/* Delete all temporary arrays and free memory */
	delete [] pq_WP_overlap_array;
	delete [] Atilde_store;
	delete [] SixJ_array;
	delete [] ClebschGordan_data;
	delete [] prefac_L_array;
	delete [] prefac_l_array;
	delete [] gsl_Plm_1_array;
	delete [] gsl_Plm_2_array;
	delete [] gsl_Plm_3_array;
	delete [] sin_phi_array;
	delete [] cos_phi_array;
  	delete [] phi_array;
  	delete [] wphi_array;

	if (use_dense_format==false){
		/* Write remaining elements to file */
		for (int thread_idx=0; thread_idx<P123_omp_num_threads; thread_idx++){
			int current_TFC = max_TFC_array[thread_idx];

			std::string thread_filename = generate_subarray_file_name(two_J_3N, P_3N,
																	  Np_WP, Nq_WP,
																	  J_2N_max,
																	  thread_idx,
																	  current_TFC,
																	  P123_folder);

			size_t num_elements_remaining = P123_dim_array_omp[thread_idx];

			/* Store array */
			store_sparse_permutation_matrix_for_3N_channel_h5(P123_val_array_omp[thread_idx],
							  								  P123_row_array_omp[thread_idx],
							  								  P123_col_array_omp[thread_idx],
							  								  num_elements_remaining,
							  								  Np_WP, p_array_WP_bounds,
							  								  Nq_WP, q_array_WP_bounds,
							  								  pw_states,
					   		  								  thread_filename,
															  false);
			/* Increment the number of files stored for current thread */
			max_TFC_array[thread_idx] += 1;

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
}

void read_and_merge_thread_files_to_single_array(double** P123_val_sparse_array,
												 int**    P123_row_array,
												 int**    P123_col_array,
												 size_t&  P123_dim,
												 int*	  max_TFC_array,
												 int	  P123_omp_num_threads,
												 int      Np_WP, double *p_array_WP_bounds,
												 int      Nq_WP, double *q_array_WP_bounds,
												 int      Nx, double* x_array, double* wx_array,
												 int      Nphi,
												 int      J_2N_max,
												 pw_3N_statespace pw_states,
												 std::string P123_folder){

	int Nalpha 	 = pw_states.Nalpha;
	int two_J_3N = pw_states.two_J_3N_array[0];
	int P_3N  	 = pw_states.P_3N_array[0];

	size_t P123_dense_dim = Nalpha * Nq_WP * Np_WP;

	std::vector<std::string> TF_filenames_vector;
	
	printf("   - Checking parallel thread files exist and reading dimensions \n"); fflush(stdout);
	/* Read dimensions from sparse subarray files */
	for (int thread_idx=0; thread_idx<P123_omp_num_threads; thread_idx++){
		int TFC_max = max_TFC_array[thread_idx];
		for (int current_TFC=0; current_TFC<TFC_max; current_TFC++){
			/* Generate filename for current thread_idx and TFC */
			std::string thread_filename = generate_subarray_file_name(two_J_3N, P_3N,
																	  Np_WP, Nq_WP,
																	  J_2N_max,
																	  thread_idx,
																	  current_TFC,
																	  P123_folder);
			
			/* Verify that file exists */
			if (std::filesystem::exists(thread_filename)==false){
				printf("     - WARNING file not found: %s \n", thread_filename.c_str()); 
				printf("       Are you certain input P123 thread files are complete? \n");
				printf("       Program will ignore file and continue ... \n");
				fflush(stdout);
				continue;
			}

			/* Append filename since it exists */
			TF_filenames_vector.push_back(thread_filename);
				
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
	for (size_t idx_TF=0; idx_TF<TF_filenames_vector.size(); idx_TF++){
		
		/* Retrieve filename from vector */
		std::string thread_filename = TF_filenames_vector[idx_TF];
				
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
													     pw_states,
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
														pw_3N_statespace pw_states,
														run_params run_parameters,
														std::string P123_folder){

	int Nalpha = pw_states.Nalpha;

	bool print_content = true;

	/* Tells program whether to run calculation, or read P123-subarrays
	 * from thread-files that failed to merge in some previous run */
	bool P123_recovery = run_parameters.P123_recovery;

	/* Array containing maximum number of sub-arrays used be each thread.
	 * Each sub-array is delt a thread file count (TFC), uniquely identifying it. */
	int* max_TFC_array = new int [run_parameters.P123_omp_num_threads];

	/* Either use a dense-block or allocate space dynamically
	 * as non-zero elements are found. The 2nd approach is a bit slower,
	 * but necessary due to the typical size of the dense-block being unhandleable. */
	bool use_dense_format = false;
	
	/* Dense and sparse array dimensions */
	size_t P123_array_dense_dim  = Nalpha * Np_WP * Nq_WP;
	size_t P123_array_dim        = 0;

	/* Temporary dense matrix to be filled in subroutine if use_dense_format=true */
	double* P123_val_dense_array = NULL;

	/* Here we initialize the calculation of non-zero P123-elements. They are
	 * calculated in chunks of 1 GB by EACH thread, which are then stored to file.
	 * This is to prevent memory overflow. The elements are then congruated in the next
	 * step. max_TFC_array is used to keep track of the subarray files. */
	if (P123_recovery==false){
		if (print_content){
			printf(" - Begin calculating 3N-channel permutation matrix ... \n");
		}
		calculate_permutation_elements_for_3N_channel(&P123_val_dense_array,
													  max_TFC_array,
													  use_dense_format,
													  production_run,
								  					  Np_WP, p_array_WP_bounds,
								  					  Nq_WP, q_array_WP_bounds,
								  					  Nx, x_array, wx_array,
													  Nphi,
													  J_2N_max,
								  					  pw_states,
													  run_parameters,
													  P123_folder);
		if (print_content){
			printf("   - Calculation finished. \n");
			fflush(stdout);
		}
	}
	else{
		for (int idx_thread=0; idx_thread<run_parameters.P123_omp_num_threads; idx_thread++){
			max_TFC_array[idx_thread] = run_parameters.max_TFC;
		}
	}
	
	/* Merge thread-subarrays of P123 to one large array and save to disk. */
	if (print_content){
		printf(" - Merging P123 parallel-distributed arrays into a single COO-format structure ... \n");
		fflush(stdout);
	}
	read_and_merge_thread_files_to_single_array(P123_sparse_val_array,
												P123_sparse_row_array,
												P123_sparse_col_array,
												P123_array_dim,
												max_TFC_array,
												run_parameters.P123_omp_num_threads,
												Np_WP, p_array_WP_bounds,
								  				Nq_WP, q_array_WP_bounds,
								  				Nx, x_array, wx_array,
												Nphi,
												J_2N_max,
								  				pw_states,
												P123_folder);
	if (print_content){
		printf("   - Merge finished. \n");
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

void fill_P123_arrays(double** P123_sparse_val_array,
					  int**    P123_sparse_row_array,
					  int**    P123_sparse_col_array,
					  size_t&  P123_sparse_dim,
					  bool     production_run,
					  fwp_statespace fwp_states,
					  pw_3N_statespace pw_states,
					  run_params run_parameters,
					  std::string P123_folder){
	
	int     Np_WP 			  = fwp_states.Np_WP;
	int     Nq_WP 			  = fwp_states.Nq_WP;
	double* p_array_WP_bounds = fwp_states.p_WP_array;
	double* q_array_WP_bounds = fwp_states.q_WP_array;

	int J_2N_max = pw_states.J_2N_max;
	int two_J_3N = pw_states.two_J_3N_array[0];
	int P_3N 	 = pw_states.P_3N_array[0];

	int Nphi = run_parameters.Nphi;
	int Nx   = run_parameters.Nx;
	
	/* Default filename for current chn_3N - used for storage and reading P123 */
	std::string P123_filename = run_parameters.P123_folder + "/" + "P123_sparse_JP_"
									+ to_string(two_J_3N) + "_" + to_string(P_3N)
									+ "_Np_" + to_string(Np_WP) + "_Nq_" + to_string(Nq_WP)
									+ "_J2max_" + to_string(J_2N_max) + ".h5";
									
	if (run_parameters.calculate_and_store_P123){
		double* x_array  = new double [Nx];
		double* wx_array = new double [Nx];
		gauss(x_array, wx_array, Nx);

		printf("Calculating P123 ... \n");
		auto timestamp_P123_calc_start = chrono::system_clock::now();
		calculate_permutation_matrices_for_all_3N_channels(P123_sparse_val_array,
														   P123_sparse_row_array,
														   P123_sparse_col_array,
														   P123_sparse_dim,
														   run_parameters.production_run,
														   Np_WP, p_array_WP_bounds,
														   Nq_WP, q_array_WP_bounds,
														   Nx, x_array, wx_array,
														   Nphi,
														   J_2N_max,
														   pw_states,
														   run_parameters,
														   run_parameters.P123_folder);
		auto timestamp_P123_calc_end = chrono::system_clock::now();
		chrono::duration<double> time_P123_calc = timestamp_P123_calc_end - timestamp_P123_calc_start;
		printf(" - Done. Time used: %.6f\n", time_P123_calc.count());

		printf("Storing P123 to h5 ... \n");
		auto timestamp_P123_store_start = chrono::system_clock::now();
		store_sparse_permutation_matrix_for_3N_channel_h5(*P123_sparse_val_array,
														  *P123_sparse_row_array,
														  *P123_sparse_col_array,
														  P123_sparse_dim,
														  Np_WP, p_array_WP_bounds,
														  Nq_WP, q_array_WP_bounds,
														  pw_states,
														  P123_filename,
														  true);
		auto timestamp_P123_store_end = chrono::system_clock::now();
		chrono::duration<double> time_P123_store = timestamp_P123_store_end - timestamp_P123_store_start;
		printf(" - Done. Time used: %.6f\n", time_P123_store.count());
	}
	else if (run_parameters.solve_faddeev){
		printf("Reading P123 from h5 ... \n");

		auto timestamp_P123_read_start = chrono::system_clock::now();
		read_sparse_permutation_matrix_for_3N_channel_h5(P123_sparse_val_array,
														 P123_sparse_row_array,
														 P123_sparse_col_array,
														 P123_sparse_dim,
														 Np_WP, p_array_WP_bounds,
														 Nq_WP, q_array_WP_bounds,
														 pw_states,
														 P123_filename,
														 true);
		auto timestamp_P123_read_end = chrono::system_clock::now();
		chrono::duration<double> time_P123_read = timestamp_P123_read_end - timestamp_P123_read_start;
		printf(" - Done. Time used: %.6f\n", time_P123_read.count());
		
		///* OLD CODE SNIPPET TO DOUBLE-CHECK ANY OPTIMIZATIONS OF P123 CALCULATON
		// * SHOULD BE MOVED TO A UNIT-TEST */
		//if (P123_sparse_dim_t==P123_sparse_dim){
		//	//int row_idx = 0;
		//	for (int idx=0; idx<P123_sparse_dim; idx++){
		//		//if (P123_sparse_row_array[idx]==row_idx){
		//		bool check1 = (abs(P123_sparse_val_array_t[idx]-P123_sparse_val_array[idx])>1e-15);
		//		bool check2 = (P123_sparse_row_array_t[idx]!=P123_sparse_row_array[idx]);
		//		bool check3 = (P123_sparse_col_array_t[idx]!=P123_sparse_col_array[idx]);
		//		if (check1||check2||check3){
		//			std::cout << "Value wrong, idx: " << idx << std::endl;
		//			std::cout << "BM val:   " << P123_sparse_val_array_t[idx] << std::endl;
		//			std::cout << "BM row:   " << P123_sparse_row_array_t[idx] << std::endl;
		//			std::cout << "BM col:   " << P123_sparse_col_array_t[idx] << std::endl;
		//			std::cout << "Prog val: " << P123_sparse_val_array[idx] << std::endl;
		//			std::cout << "Prog row: " << P123_sparse_row_array[idx] << std::endl;
		//			std::cout << "Prog col: " << P123_sparse_col_array[idx] << std::endl;
		//			raise_error("element mismatch");
		//		}
		//		//}
		//	}
		//}
		//else{
		//	std::cout << "BM dim:   " << P123_sparse_dim_t << std::endl;
		//	std::cout << "Prog dim: " << P123_sparse_dim << std::endl;
		//	raise_error("dim not right");
		//}
	}
}
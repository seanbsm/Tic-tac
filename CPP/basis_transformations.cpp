
#include "basis_transformations.h"

/* This function calculates all CG coefficients for a given j1 and j2,
 * where j1 is half-integer and j2 is an integer */
void calculate_CG_array_for_half_integers(double* CG_array,
										  int two_j1_max,
										  int two_j2_max){

	int max_num_j1_vals = (two_j1_max+1)/2;
	int max_num_j2_vals =  two_j2_max/2;
	int max_num_j3_vals = (two_j2_max + two_j1_max - abs(two_j2_max-two_j1_max) + 2)/2;

	double CG_coeff_sum = 0;
	/* two_j1 - loop (note starting value 1 and step-length 2 since j1 is a half-integer) */
	for (int two_j1=1; two_j1<two_j1_max+1; two_j1+=2){
		/* two_j2 - loop (note starting value 0 and step-length 2 since j2 is an integer) */
		for (int two_j2=0; two_j2<two_j2_max+1; two_j2+=2){
		
			int two_j3_min = (int) abs(two_j1 - two_j2);
			int two_j3_max = two_j1 + two_j2;
			/* two_j3 - loop (note starting value 1 and step-length 2 since j1 is a half-integer) */
			for (int two_j3=two_j3_min; two_j3<two_j3_max+1; two_j3+=2){

				CG_coeff_sum = 0;
				/* two_m_j1 - loop */
				for (int two_m_j1=-two_j1; two_m_j1<two_j1+1; two_m_j1+=2){
					/* two_m_j2 - loop */
					for (int two_m_j2=-two_j2; two_m_j2<two_j2+1; two_m_j2+=2){
						/* two_m_j3 = two_m_j2 + two_m_j1 (CG-coeff is zero otherwise)*/ 
						int two_m_j3 = two_m_j2 + two_m_j1;

						CG_coeff_sum += clebsch_gordan(two_j1, two_j2, two_j3, two_m_j1, two_m_j2, two_m_j3);
					}
				}

				/* Note the division by two - necessary when using quantum numbers for index-value calculation */
				int CG_array_idx = ((two_j1-1)*two_j1_max*two_j2_max + two_j2*two_j3_max + two_j3-1) / 2;

				/* Set sum in array */
				CG_array[CG_array_idx] = CG_coeff_sum;
			}
		}
	}
}

/* This function calculates all CG coefficients for a given j1 and j2,
 * where both j1 and j2 are whole integers */
void calculate_CG_array_for_whole_integers(double* CG_array,
										   int j1_max,
										   int j2_max){

	int max_num_j1_vals = j1_max;
	int max_num_j2_vals = j2_max;
	int max_num_j3_vals = j2_max + j1_max - abs(j2_max-j1_max) + 1;

	double CG_coeff_sum = 0;
	/* j1 - loop*/
	for (int j1=0; j1<j1_max+1; j1++){
		/* j2 - loop */
		for (int j2=0; j2<j2_max+1; j2++){
		
			int j3_min = (int) abs(j1 - j2);
			int j3_max = j1 + j2;
			/* j3 - loop */
			for (int j3=0; j3<j3_max+1; j3++){

				CG_coeff_sum = 0;
				/* m_j1 - loop */
				for (int m_j1=-j1; m_j1<j1+1; m_j1++){
					/* m_j2 - loop */
					for (int m_j2=-j2; m_j2<j2+1; m_j2++){
						/* m_j3 = m_j2 + m_j1 (CG-coeff is zero otherwise)*/ 
						int m_j3 = m_j2 + m_j1;

						CG_coeff_sum += clebsch_gordan(2*j1, 2*j2, 2*j3, 2*m_j1, 2*m_j2, 2*m_j3);
					}
				}

				/* Note the division by two - necessary when using quantum numbers for index-value calculation */
				int CG_array_idx = j1*j1_max*j2_max + j2*j3_max + j3;

				/* Set sum in array */
				CG_array[CG_array_idx] = CG_coeff_sum;
			}
		}
	}
}

void precalculate_clebsch_gordan_coefficients(double* CG_coeffs_array_J,
											  double* CG_coeffs_array_L,
											  double* CG_coeffs_array_S,
											  double* CG_coeffs_array_T,
											  int two_J_1N_max,
											  int J_2N_max){

	/* Total angular momentum J - couplings */
	int max_num_J_1N_vals = (two_J_1N_max+1)/2;
	int max_num_J_2N_vals = J_2N_max;
	int max_num_J_3N_vals = (2*J_2N_max+two_J_1N_max - abs(2*J_2N_max-two_J_1N_max) + 2)/2;

	double* CG_J_array = new double [max_num_J_1N_vals * max_num_J_2N_vals * max_num_J_3N_vals];

	calculate_CG_array_for_half_integers(CG_J_array,
										 two_J_1N_max,
										 2*J_2N_max);
	
	/* Total spin S - couplings */
	int max_num_S_1N_vals = 1;  // S_1N = 1/2
	int max_num_S_2N_vals = 2;  // S_2N = 0 or 1
	int max_num_S_3N_vals = 2;  // S_3N = 1/2 or 3/2

	double* CG_S_array = new double [max_num_S_1N_vals * max_num_S_2N_vals * max_num_S_3N_vals];

	calculate_CG_array_for_half_integers(CG_S_array,
										 1,
										 2);

	/* Total spin T - couplings (identical to the S-calculation above) */
	int max_num_T_1N_vals = 1;  // T_1N = 1/2
	int max_num_T_2N_vals = 2;  // T_2N = 0 or 1
	int max_num_T_3N_vals = 2;  // T_3N = 1/2 or 3/2

	double* CG_T_array = new double [max_num_T_1N_vals * max_num_T_2N_vals * max_num_T_3N_vals];

	calculate_CG_array_for_half_integers(CG_T_array,
										 1,
										 2);
	
	/* Orbital angular momentum - couplings (distinct from the calculations above as j1,j2,j3 are all whole integers) */
	int max_num_L_1N_vals = max_num_J_1N_vals+1;
	int max_num_L_2N_vals = max_num_J_2N_vals+1;
	int max_num_L_3N_vals = max_num_J_3N_vals+2;

	double* CG_L_array = new double [max_num_L_1N_vals * max_num_L_2N_vals * max_num_L_3N_vals];

	calculate_CG_array_for_whole_integers(CG_T_array,
										  (two_J_1N_max+1)/2,
										  J_2N_max+1);
}

int calculate_SH_array_idx(int p, int l, int m, int l_max){
	/* p = 0, 1, ..., num_angles-1
	 * l = 0, 1, ..., l_max
	 * m = -l, -l+1, ... 0, ..., l-1, l */
	
	return p*(l_max+1)*(l_max+1) + l*l + m + l;
}

void generate_spherical_harmonics_array(double **SH_array, double* angles_array, int l_max, int num_angles){
	
	int array_size_per_angle = gsl_sf_legendre_array_n(l_max);
	
	*SH_array = new double [num_angles*(l_max+1)*(l_max+1)];

	double gsl_SH_array [array_size_per_angle];
	
	for(int p=0; p<num_angles; p++){
		double c = cos( angles_array[p]*deg_to_rad );
		
		/* The 4th argument is set to -1, which means we use the Condon-Shortley phase factor (-1)^m */
		gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM, l_max, c,-1, gsl_SH_array);
		
		for (int l=0; l<l_max+1; l++){
			for (int m=0; m<l+1; m++){
				
				int idx_gsl = gsl_sf_legendre_array_index(l, m);
				
				/* Positive m */
				int idx_SH  = calculate_SH_array_idx(p, l, m, l_max);
				*SH_array[idx_SH] = (double)gsl_SH_array[idx_gsl];
				
				/* Negative m */
				if (m!=0){
					double factor = tgamma(l-m+1)/tgamma(l+m+1);
					factor *= 1 - 2*(m%2);
					
					int idx_SH  = calculate_SH_array_idx(p, l, -m, l_max);
					*SH_array[idx_SH] = factor * (floatType)gsl_SH_array[idx_gsl];
				}
			}
		}
	}
}

//cdouble calculate_U_channel_spin_element(cdouble* U_array,
//                                         int two_J_3N, int two_T_3N, int J_2N, int S_2N, int T_2N, int L_2N,
//                                         int alpha_idx_r, int two_Sigma_r,
//                                         int alpha_idx_c, int two_Sigma_c,
//                                         int Nalpha, int Np_WP, int Nq_WP,
//                                         int* L_2N_array,
//                                         int* S_2N_array,
//                                         int* J_2N_array,
//                                         int* T_2N_array,
//                                         int* L_1N_array,
//                                         int* two_J_1N_array,
//                                         int* two_J_3N_array,
//                                         int* two_T_3N_array){
//    
//    int U_dim = Nalpha * Np_WP * Nq_WP;
//
//    int two_J_1N_r = two_J_1N_array[alpha_idx_r];
//    int two_J_1N_c = two_J_1N_array[alpha_idx_c];
//
//    int L_1N_r = L_1N_array[alpha_idx_r];
//    int L_1N_c = L_1N_array[alpha_idx_c];
//
//    /* This is 1 if it's odd, 0 otherwise */
//    bool odd = ((two_J_3N-two_J_1N_r)/2 )%2;
//    double coupling_coeffs = 1;
//
//    coupling_coeffs *= sqrt( (two_J_1N_r+1) * (two_Sigma_r+1) );
//    coupling_coeffs *= 1 - 2*odd;
//    coupling_coeffs *= wigner_6j(2*L_1N_r, 1, two_J_1N_r, 2*J_2N, two_J_3N, two_Sigma_r);
//
//    odd = ((two_J_3N-two_J_1N_c)/2 )%2;
//
//    coupling_coeffs *= sqrt( (two_J_1N_c+1) * (two_Sigma_c+1) );
//    coupling_coeffs *= 1 - 2*odd;
//    coupling_coeffs *= wigner_6j(2*L_1N_c, 1, two_J_1N_c, 2*J_2N, two_J_3N, two_Sigma_c);
//
//    int U_idx = (int) ( alpha_idx_c*Nq_WP*Np_WP +  q_idx_WP*Np_WP +  p_idx_WP)*U_dim +
//                        alpha_idx_r*Nq_WP*Np_WP + qp_idx_WP*Np_WP + pp_idx_WP;
//    
//    return coupling_coeffs * U_array[U_idx];
//}
//
//cdouble calculate_M_channel_spin_element(cdouble* U_array,
//                                         int two_J_3N, int two_T_3N, int J_2N, int S_2N, int T_2N, int L_2N,
//                                         int two_Sigma_r, int two_m_Sigma_r,
//                                         int two_Sigma_c, int two_m_Sigma_c,
//                                         int idx_angle,
//                                         double* SH_array,
//                                         int Nalpha,
//                                         int* L_2N_array,
//                                         int* S_2N_array,
//                                         int* J_2N_array,
//                                         int* T_2N_array,
//                                         int* L_1N_array,
//                                         int* two_J_1N_array,
//                                         int* two_J_3N_array,
//                                         int* two_T_3N_array,
//                                         int  two_J_3N_max,
//                                         int  parity_3N,
//                                         double k){
//    
//    /* Largest angular momentum possible for orbiting nucleon,
//     * obtained when J_2N=0 -> J_3N_max = (J_2N=0) + J_1N_max = L_1N_max + 1/2. 
//     * This variable is used in the spherical harmonics */
//    int l_max = (two_J_3N_max-1)/2;
//
//    /* M-element to return at end of function */
//    cdouble M_element = {0, 1.};
//    M_element        *= sqrt(M_PI)/k;
//    
//    /* Loop summation-variable */
//    cdouble sum = {0., 0.};
//    /* Loop over row-states */
//    for (int alpha_idx_r=0; alpha_idx_r<Nalpha; alpha_idx_r++){
//
//        /* Verify that this state fits the experimental setup (e.g. a deuteron target) */
//        int L_2N_r = L_2N_array[alpha_idx_r];
//        int J_2N_r = J_2N_array[alpha_idx_r];
//        int S_2N_r = S_2N_array[alpha_idx_r];
//        int T_2N_r = T_2N_array[alpha_idx_r];
//        if (J_2N_r!=J_2N || S_2N_r!=S_2N || T_2N_r!=T_2N){
//            continue;
//        }
//
//        /* Loop over column-states */
//        for (int alpha_idx_c=0; alpha_idx_c<Nalpha; alpha_idx_c++){
//
//            /* Verify that this state fits the experimental setup (e.g. a deuteron target) */
//            int L_2N_c = L_2N_array[alpha_idx_c];
//            int J_2N_c = J_2N_array[alpha_idx_c];
//            int S_2N_c = S_2N_array[alpha_idx_c];
//            int T_2N_c = T_2N_array[alpha_idx_c];
//            if (J_2N_c!=J_2N || S_2N_c!=S_2N || T_2N_c!=T_2N){
//                continue;
//            }
//
//            /* Verify the two states can couple */
//            if (abs(L_2N_r - L_2N_c)>2){
//                continue;
//            }
//
//            cdouble factors = {1., 1.};
//
//            const cdouble imag_unit = {0., 1.};
//
//            int L_1N_r = L_1N_array[alpha_idx_r];
//            int L_1N_c = L_1N_array[alpha_idx_c];
//
//            int m_Sigma_diff = (two_m_Sigma_c - two_m_Sigma_r) / 2;
//            if (abs(m_Sigma_diff) > L_1N_r){
//                continue;
//            }
//
//            /* i^(lambda' - lambda) * sqrt(2*lambda+1)*/
//            factors   *= std::pow(imag_unit, L_1N_r- L_1N_c) * sqrt(2*L_1N_c + 1);
//
//            /* Clebsch-Gordan couplings */
//            factors   *= clebsch_gordan(two_Sigma_c, 2*L_1N_c, two_J_3N, two_m_Sigma_c, 0, two_m_Sigma_c);
//            factors   *= clebsch_gordan(two_Sigma_r, 2*L_1N_r, two_J_3N, two_m_Sigma_r, 2*m_Sigma_diff, two_m_Sigma_c);
//            
//            /* Find index of pre-calculated spherical harmonics (SH) */
//            int SH_idx = calculate_SH_array_idx(idx_angle, L_1N_r, m_Sigma_diff, l_max);
//            SH_coeff   = SH_array[SH_idx];
//            factors   *= SH_coeff;
//
//            double delta_function = (two_Sigma_r == two_Sigma_c) * (L_1N_r == L_1N_c);
//
//            cdouble U_element = calculate_U_channel_spin_element(U_array,
//                                                                 two_J_3N, two_T_3N, J_2N, S_2N, T_2N, L_2N,
//                                                                 L_1N_r, two_Sigma_r,
//                                                                 L_1N_c, two_Sigma_c,
//                                                                 Nalpha,
//                                                                 L_2N_array,
//                                                                 S_2N_array,
//                                                                 J_2N_array,
//                                                                 T_2N_array,
//                                                                 L_1N_array,
//                                                                 two_J_1N_array,
//                                                                 two_J_3N_array,
//                                                                 two_T_3N_array);
//
//            sum += factors*(delta_function - U_element);
//
//            /* The 4th argument is set to -1, which means we use the Condon-Shortley phase factor (-1)^m */
//		    //gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM, lmax, c,-1, spherical_harmonis_array);
//        }
//    }
//
//    M_element *= sum;
//
//    return M_element;
//}
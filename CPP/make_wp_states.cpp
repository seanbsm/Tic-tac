#include "make_wp_states.h"

double p_normalization(double p0, double p1){
	//double Ep0 = p0*p0/MN;
	//double Ep1 = p1*p1/MN;
	//return sqrt(Ep1-Ep0);	// energy WPs
	return sqrt(p1-p0);	// momentum WPs
}
double q_normalization(double q0, double q1){
	//double Eq0 = com_q_momentum_to_com_energy(q0);
	//double Eq1 = com_q_momentum_to_com_energy(q1);
	//return sqrt(Eq1-Eq0);	// energy WPs
	//return sqrt( (q1*q1-q0*q0)/MN );	// energy WPs (old)
	return sqrt(q1-q0);	// momentum WPs
}

double p_weight_function(double p){
	//return sqrt(2*p/MN);	// energy WPs
	return 1;				// momentum WPs
}
double q_weight_function(double q){
	//return sqrt(2*q/MN);	// energy WPs
	return 1;				// momentum WPs
}

void make_chebyshev_distribution(int N_WP,
								 double* boundary_array,
								 double scale,
								 int 	sparseness_degree){
	double tan_term = 0;
	double boundary = 0;
	
	for (int i=1; i<N_WP+1; i++){
		tan_term = tan( (2*i-1)*M_PI/(4*N_WP) ); 

		boundary = scale*pow(tan_term, sparseness_degree);
		
		/* Momentum distribution */
		boundary_array[i] = boundary;
		/* Energy distribution */
		//boundary_array[i] = com_energy_to_com_q_momentum(boundary);
	}
	boundary_array[0] = 0.0;
}

void make_p_bin_grid(int Np_WP, double* p_WP_array, run_params run_parameters){
	//double min_boundary 	 = 1e-5;
	//double max_boundary 	 = 100;
	//double scale = 0.5*(max_boundary - min_boundary);
	double scale = run_parameters.chebyshev_s;
	int    sparseness_degree = run_parameters.chebyshev_t;

	make_chebyshev_distribution(Np_WP, p_WP_array,
								scale,
								sparseness_degree);
}
void make_q_bin_grid(int Nq_WP, double* q_WP_array, run_params run_parameters){
	//double min_boundary 	 = 1e-5;
	//double max_boundary 	 = 100;
	//double scale = 0.5*(max_boundary - min_boundary);
	double scale = run_parameters.chebyshev_s;
	int    sparseness_degree = run_parameters.chebyshev_t;

	make_chebyshev_distribution(Nq_WP, q_WP_array,
								scale,
								sparseness_degree);
}

void make_p_bin_quadrature_grids(int Np_WP, double* p_WP_array,
								 int Np_per_WP, double* p_array, double* wp_array){
	for (int idx_bin=0; idx_bin<Np_WP; idx_bin++){
		double bin_lower_bound = p_WP_array[idx_bin];
		double bin_upper_bound = p_WP_array[idx_bin+1];

		double*  p_array_ptr =  &p_array[idx_bin*Np_per_WP];
		double* wp_array_ptr = &wp_array[idx_bin*Np_per_WP];
		for (int idx_p=0; idx_p<Np_per_WP; idx_p++){
			gauss(p_array_ptr, wp_array_ptr, Np_per_WP);
			updateRange_a_b(p_array_ptr, wp_array_ptr, bin_lower_bound, bin_upper_bound, Np_per_WP);
		}
	}
}
void make_q_bin_quadrature_grids(int Nq_WP, double* q_WP_array,
								 int Nq_per_WP, double* q_array, double* wq_array){
	for (int idx_bin=0; idx_bin<Nq_WP; idx_bin++){
		double bin_lower_bound = q_WP_array[idx_bin];
		double bin_upper_bound = q_WP_array[idx_bin+1];

		double*  q_array_ptr =  &q_array[idx_bin*Nq_per_WP];
		double* wq_array_ptr = &wq_array[idx_bin*Nq_per_WP];
		for (int idx_p=0; idx_p<Nq_per_WP; idx_p++){
			gauss(q_array_ptr, wq_array_ptr, Nq_per_WP);
			updateRange_a_b(q_array_ptr, wq_array_ptr, bin_lower_bound, bin_upper_bound, Nq_per_WP);
		}
	}
}
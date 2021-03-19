#include "make_wp_states.h"

double p_normalization(double p0, double p1){
	//return sqrt( (p1*p1-p0*p0)/MN );
	return sqrt(p1-p0);
}
double q_normalization(double q0, double q1){
	//return sqrt( (q1*q1-q0*q0)/MN );
	return sqrt(q1-q0);
}

double p_weight_function(double p){
	//return sqrt(p/MN);
	return 1;
}
double q_weight_function(double q){
	//return sqrt(q/MN);
	return 1;
}

void make_chebyshev_distribution(int N_WP,
								 double* boundary_array,
								 double min_boundary,
								 double max_boundary,
								 int 	sparseness_degree){
	double tan_term;
	
	double scale = 0.5*(max_boundary - min_boundary);
	
	for (int i=0; i<N_WP+1; i++){
		tan_term = tan( (2*i)*M_PI/(4*(N_WP+1)) );
		
		//boundary_array[i] = i;
		boundary_array[i] = scale*pow(tan_term, sparseness_degree);
	}
	//boundary_array[0] = 0.001;
}

void make_p_bin_grid(int Np_WP, double* p_WP_array){
	double min_boundary 	 = 1e-5;
	double max_boundary 	 = 100;
	int    sparseness_degree = 2;

	make_chebyshev_distribution(Np_WP, p_WP_array,
								min_boundary, max_boundary,
								sparseness_degree);
}
void make_q_bin_grid(int Nq_WP, double* q_WP_array){
	double min_boundary      = 1e-5;
	double max_boundary      = 100;
	int    sparseness_degree = 2;

	make_chebyshev_distribution(Nq_WP, q_WP_array,
								min_boundary, max_boundary,
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
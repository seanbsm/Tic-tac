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
								 double	sparseness_degree){
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

void make_p_bin_grid(fwp_statespace& fwp_states, run_params run_parameters){
	if (run_parameters.p_grid_type=="chebyshev"){
		double scale 			 = run_parameters.chebyshev_s;
		double sparseness_degree = run_parameters.chebyshev_t;
		make_chebyshev_distribution(fwp_states.Np_WP, fwp_states.p_WP_array,
									scale,
									sparseness_degree);
	}
	else if (run_parameters.p_grid_type=="custom"){
		read_WP_boundaries_from_txt(fwp_states.p_WP_array, fwp_states.Np_WP, run_parameters.p_grid_filename);
	}
	else{
		raise_error("Unknown p-momentum gridtype specified.");
	}
}
void make_q_bin_grid(fwp_statespace& fwp_states, run_params run_parameters){
	if (run_parameters.q_grid_type=="chebyshev"){
		double scale 			 = run_parameters.chebyshev_s;
		double sparseness_degree = run_parameters.chebyshev_t;
		make_chebyshev_distribution(fwp_states.Nq_WP, fwp_states.q_WP_array,
									scale,
									sparseness_degree);
	}
	else if (run_parameters.p_grid_type=="custom"){
		read_WP_boundaries_from_txt(fwp_states.q_WP_array, fwp_states.Nq_WP, run_parameters.q_grid_filename);
	}
	else{
		raise_error("Unknown q-momentum gridtype specified.");
	}
}

void make_p_bin_quadrature_grids(fwp_statespace& fwp_states){
	int 	Np_WP		= fwp_states.Np_WP;
	double* p_WP_array	= fwp_states.p_WP_array;
	int 	Np_per_WP	= fwp_states.Np_per_WP;
	double* p_array		= fwp_states.p_array;
	double* wp_array	= fwp_states.wp_array;

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
void make_q_bin_quadrature_grids(fwp_statespace& fwp_states){
	int 	Nq_WP		= fwp_states.Nq_WP;
	double* q_WP_array	= fwp_states.q_WP_array;
	int 	Nq_per_WP	= fwp_states.Nq_per_WP;
	double* q_array		= fwp_states.q_array;
	double* wq_array	= fwp_states.wq_array;

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

void make_fwp_statespace(fwp_statespace& fwp_states, run_params run_parameters){
	printf("Constructing wave-packet (WP) state space ... \n");

	/* Copy space dimensions from input struct */
	fwp_states.Np_WP 	= run_parameters.Np_WP;
	fwp_states.Nq_WP 	= run_parameters.Nq_WP;
	fwp_states.Np_per_WP = run_parameters.Np_per_WP;
	fwp_states.Nq_per_WP = run_parameters.Nq_per_WP;

	/* Make bin boundaries */
	printf(" - Constructing wave-packet (WP) p-momentum bin boundaries ... \n");
	fwp_states.p_WP_array = new double [fwp_states.Np_WP+1];
	make_p_bin_grid(fwp_states, run_parameters);
	printf("   - Done \n");
	printf(" - Constructing wave-packet (WP) q-momentum bin boundaries ... \n");
	fwp_states.q_WP_array = new double [fwp_states.Nq_WP+1];
	make_q_bin_grid(fwp_states, run_parameters);
	printf("   - Done \n");

	/* Make Gauss-Legendre quadrature meshes inside each bin */
	printf(" - Constructing p quadrature mesh per WP, for all WPs ... \n");
	fwp_states.p_array  = new double [fwp_states.Np_per_WP * fwp_states.Np_WP];
	fwp_states.wp_array = new double [fwp_states.Np_per_WP * fwp_states.Np_WP];
	make_p_bin_quadrature_grids(fwp_states);
	printf("   - Done \n");
	printf(" - Constructing q quadrature mesh per WP, for all WPs ... \n");
	fwp_states.q_array  = new double [fwp_states.Nq_per_WP * fwp_states.Nq_WP];
	fwp_states.wq_array = new double [fwp_states.Nq_per_WP * fwp_states.Nq_WP];
	make_q_bin_quadrature_grids(fwp_states);
	printf("   - Done \n");

	/* Store boundaries for post-processing of output  */
	printf(" - Storing q boundaries to CSV-file ... \n");
	std::string q_boundaries_filename = run_parameters.output_folder + "/" + "q_boundaries_Nq_" + std::to_string(fwp_states.Nq_WP) + ".csv";
	store_q_WP_boundaries_csv(fwp_states, q_boundaries_filename);
	printf("   - Done \n");

	printf(" - Done \n");
}
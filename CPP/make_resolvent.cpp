
#include "make_resolvent.h"

double heaviside_step_function(double val){
	if (val<0){
		return 0;
	}
	else{
		return 1;
	}
}

/* See header-file, commentary (A), for explanation of notation and equations */
cdouble resolvent_bound_continuum(double E, double Eb,
								  double q_bin_upper,
								  double q_bin_lower){

	double mu1 = Mn*(Mn+Mp+Eb)/(Mn + Mn + Mp + Eb);

	/* Bin boundaries in energy */
	double Eq_lower = q_bin_lower*q_bin_lower/(2*mu1);
	double Eq_upper = q_bin_upper*q_bin_upper/(2*mu1);

	/* Energy width of SWPs (D is short for Delta) for q momenta */
	double Dq = Eq_upper - Eq_lower;

	/* Real part of the BC resolvent */
	double Re_R = std::log( std::abs( (Eq_lower + Eb - E)/(Eq_upper + Eb - E) ) ) / (Dq); 
	
	/* Imaginary part of the BC resolvent */
	double Im_R = (  heaviside_step_function( Eq_upper + Eb - E )
				   - heaviside_step_function( Eq_lower + Eb - E ) ) * (-M_PI)/ (Dq); 
	
	/* Return the complex BC resolvent term Q */
	return {Re_R, Im_R};
}

/* See header-file, commentary (A), for explanation of notation and equations */
cdouble resolvent_continuum_continuum(double E, double Eb,
									  double q_bin_upper,
									  double q_bin_lower,
									  double e_bin_upper,
									  double e_bin_lower){

	double mu1 = Mn*(Mn+Mp+Eb)/(Mn + Mn + Mp + Eb);

	/* Bin boundaries in energy */
	double Eq_lower = q_bin_lower*q_bin_lower/(2*mu1);
	double Eq_upper = q_bin_upper*q_bin_upper/(2*mu1);

	double Ep_lower = e_bin_lower;
	double Ep_upper = e_bin_upper;

	double Ep = 0.5*(Ep_lower + Ep_upper);
	double Eq = 0.5*(Eq_lower + Eq_upper);

	/* Energy widths of SWPs (D is short for Delta) for p and q momenta */
	double Dq = Eq_upper - Eq_lower;
	double Dp = Ep_upper - Ep_lower;

	/* Temporary variables */
	double DM = 0.5*(Dp - Dq);
	double DP = 0.5*(Dp + Dq);
	double D  = Ep + Eq - E;
	
	/* Real part of the CC resolvent */
	double Re_Q = (  (D+DM) * std::log( std::abs(D+DM) )
				   + (D-DM) * std::log( std::abs(D-DM) )
				   - (D+DP) * std::log( std::abs(D+DP) )
				   - (D-DP) * std::log( std::abs(D-DP) ) ) / (Dp*Dq); 
	
	/* Imaginary part of the CC resolvent */
	double Im_Q = (  (D+DM) * heaviside_step_function( D+DM )
				   + (D-DM) * heaviside_step_function( D-DM )
				   - (D+DP) * heaviside_step_function( D+DP )
				   - (D-DP) * heaviside_step_function( D-DP ) ) * M_PI/ (Dp*Dq); 
	
	/* Return the complex CC resolvent term Q */
	return {Re_Q, Im_Q};
}

void calculate_resolvent_array_in_SWP_basis(cdouble* G_array,
											double   E,
											swp_statespace swp_states,
											pw_3N_statespace pw_states,
											run_params run_parameters){
	
	int 	Np_WP			 = swp_states.Np_WP;
	int     Nq_WP			 = swp_states.Nq_WP;
	double* e_SWP_unco_array = swp_states.e_SWP_unco_array;
	double* e_SWP_coup_array = swp_states.e_SWP_coup_array;
	double* q_WP_array		 = swp_states.q_WP_array;

	int  Nalpha			= pw_states.Nalpha;
	int* L_2N_array		= pw_states.L_2N_array;
	int* S_2N_array		= pw_states.S_2N_array;
	int* J_2N_array		= pw_states.J_2N_array;
	int* T_2N_array		= pw_states.T_2N_array;
	int* two_T_3N_array = pw_states.two_T_3N_array;

	/* This test will be reused several times */
	bool tensor_force_true = (run_parameters.tensor_force==true);

	/* Pointer to either p_SWP_unco_array or p_SWP_coup_array,
	 * which is determined by whether the channel is coupled or not */
	double* e_SWP_array_ptr = NULL;

	/* Loop over states along resolvent diagonal */
	for (int idx_alpha=0; idx_alpha<Nalpha; idx_alpha++){
	
		int L = L_2N_array[idx_alpha];
		int S = S_2N_array[idx_alpha];
		int J = J_2N_array[idx_alpha];
		int T = T_2N_array[idx_alpha];

		int two_T_3N = two_T_3N_array[idx_alpha];

		/* Detemine if this is a coupled channel.
		 * !!! With isospin symmetry-breaking we count 1S0 as a coupled matrix via T_3N-coupling !!! */
		bool coupled_channel = false;
		bool state_1S0 = (S==0 && J==0 && L==0);
		bool coupled_via_L_2N = (tensor_force_true && L!=J && J!=0);
		bool coupled_via_T_3N = (state_1S0==true && run_parameters.isospin_breaking_1S0==true);
		if (coupled_via_L_2N && coupled_via_T_3N){
			raise_error("Warning! Code has not been written to handle isospin-breaking in coupled channels!");
		}
		if (coupled_via_L_2N || coupled_via_T_3N){ // This counts 3P0 as uncoupled; used in matrix structure
			coupled_channel  = true;
		}

		if (coupled_channel){
			int chn_2N_idx = unique_2N_idx(L, S, J, T, coupled_channel, run_parameters);
			if (coupled_via_L_2N){
				if (L<J){
					e_SWP_array_ptr = &e_SWP_coup_array[chn_2N_idx * 2*(Np_WP+1)];
				}
				else{
					e_SWP_array_ptr = &e_SWP_coup_array[chn_2N_idx * 2*(Np_WP+1) + Np_WP+1];
				}
			}
			else if (coupled_via_T_3N){
				if (two_T_3N==1){
					e_SWP_array_ptr = &e_SWP_coup_array[chn_2N_idx * 2*(Np_WP+1)];
				}
				else{
					e_SWP_array_ptr = &e_SWP_coup_array[chn_2N_idx * 2*(Np_WP+1) + Np_WP+1];
				}
			}
			else{
				raise_error("Unknown coupling encountered in resolvent-calculation!");
			}
		}
		else{
			int chn_2N_idx = unique_2N_idx(L, S, J, T, coupled_channel, run_parameters);
			e_SWP_array_ptr = &e_SWP_unco_array[chn_2N_idx * (Np_WP+1)];
		}

		//printf("%d %d %d %d \n", L, S, J, T);

		/* p-momentum index loop */
		for (int idx_p_bin=0; idx_p_bin<Np_WP; idx_p_bin++){

			/* Upper and lower boundaries of current p-bin (expressed in energy) */
			double e_bin_lower = e_SWP_array_ptr[idx_p_bin    ];
			double e_bin_upper = e_SWP_array_ptr[idx_p_bin + 1];

			//printf("%.4f\n", e_bin_lower);

			/* Bound state check, given by p_bin_lower if bound state exists */
			bool bound_state_exists = false;
			double Eb = 0;
			if (e_bin_lower<0){
				Eb = e_bin_lower;
				bound_state_exists = true;
			}

			/* q-momentum index loop */            
			for (int idx_q_bin=0; idx_q_bin<Nq_WP; idx_q_bin++){

				/* Upper and lower boundaries of current q-bin */
				double q_bin_lower = q_WP_array[idx_q_bin    ];
				double q_bin_upper = q_WP_array[idx_q_bin + 1];

				cdouble R = {0, 0};
				cdouble Q = {0, 0};
				if (bound_state_exists){    // Calculate bound-continuum (BC) resolvent part R
					R = resolvent_bound_continuum(E, Eb,
												  q_bin_upper,
												  q_bin_lower);
				}
				else{                       // Calculate continuum-continuum (CC) resolvent part Q
					Q = resolvent_continuum_continuum(E, Eb,
													  q_bin_upper,
													  q_bin_lower,
													  e_bin_upper,
													  e_bin_lower);
					//std::cout << Q << std::endl;
				}

				/* Use identical indexing as used in permutation matrix */
				int G_idx = idx_alpha*Nq_WP*Np_WP + idx_q_bin*Np_WP + idx_p_bin;
				G_array[G_idx] = R + Q;

				//if (e_bin_lower<0){
				//	std::cout << R << " " << Q << std::endl;
				//}
				//std::cout << R << " " << Q << std::endl;
			}
		}
	}
}
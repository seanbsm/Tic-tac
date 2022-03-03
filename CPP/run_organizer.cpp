
#include "run_organizer.h"

void find_on_shell_bins(solution_configuration& solve_config,
                        swp_statespace swp_states,
                        run_params run_parameters){
    
    /* Make local pointers & variables for solution configuration */
	size_t  num_T_lab	    = 0;
    /* Make local pointers & variables for swp-statespace */
    int     Nq_WP      = swp_states.Nq_WP;
    double  E_bound    = swp_states.E_bound;
    double* q_WP_array = swp_states.q_WP_array;
    
    /* Use q-momentum bin ENERGY mid-points as on-shell energies if no default input is given */
	std::vector<size_t> q_WP_idx_vec;

	/* Special condition to reduce number of on-shell calculations */
	double Eq_lower = 0;
	double Eq_upper = 0;
	double E_com    = 0;
	double q_m		= 0;
	std::vector<bool>   midpoint_idx_vector   (Nq_WP-1, false);
	std::vector<double> T_lab_midpoint_vector (Nq_WP, false);
	for (size_t q_WP_idx=0; q_WP_idx<Nq_WP; q_WP_idx++){
		Eq_lower = 0.5*(q_WP_array[q_WP_idx]   * q_WP_array[q_WP_idx])  /mu1(E_bound);
		Eq_upper = 0.5*(q_WP_array[q_WP_idx+1] * q_WP_array[q_WP_idx+1])/mu1(E_bound);
		E_com	 = 0.5*(Eq_upper + Eq_lower);
		q_m      = com_energy_to_com_q_momentum(E_com);
		T_lab_midpoint_vector[q_WP_idx] = com_momentum_to_lab_energy(q_m, E_bound);
	}

	double* energy_input_array = NULL;
	int		num_energy_input   = 0;
	read_input_energies(energy_input_array, num_energy_input, run_parameters.energy_input_file);

	//std::vector<double> T_lab_input_list = {1,2,3,4,5,9,10,13,22.7,35,53};
	//std::vector<double> T_lab_input_list = {  3,   4,   5,       6,
	//										  9,  10,  11,      12,
	//										 13,  16,  22.7,    28,
	//										 30,  35,  42,      47.5,
	//										 50,  53,  65,      93.5,
	//										146, 155, 180, 220, 240};
			
	for (size_t q_WP_idx=0; q_WP_idx<Nq_WP-1; q_WP_idx++){
		double T_lab_lower = T_lab_midpoint_vector[q_WP_idx];
		double T_lab_upper = T_lab_midpoint_vector[q_WP_idx+1];
		for (int i=0; i<num_energy_input; i++){
			double T_lab_input = energy_input_array[i];
			/* See if input energy lies between two bin mid-points */
			if (T_lab_lower<=T_lab_input && T_lab_input<=T_lab_upper){
				midpoint_idx_vector[q_WP_idx] = true;
			}
		}
	}
	/* Use on-shell midpoints to set on-shell bins */
	std::vector<bool>   bin_idx_vector   (Nq_WP, false);
	for (size_t q_WP_idx=0; q_WP_idx<Nq_WP-1; q_WP_idx++){
		if (midpoint_idx_vector[q_WP_idx]==true){
			bin_idx_vector[q_WP_idx]   = true;
			bin_idx_vector[q_WP_idx+1] = true;
		}
	}
	/* Append on-shell bin indices to q_WP_idx_vec */
	for (size_t q_WP_idx=0; q_WP_idx<Nq_WP; q_WP_idx++){
		if (bin_idx_vector[q_WP_idx]==true){
			q_WP_idx_vec.push_back(q_WP_idx);
		}
	}
	num_T_lab = q_WP_idx_vec.size();

	solve_config.T_lab_array = new double [num_T_lab];

	for (size_t Tlab_idx=0; Tlab_idx<num_T_lab; Tlab_idx++){
		size_t q_WP_idx = q_WP_idx_vec[Tlab_idx];
		double Eq_lower = 0.5*(q_WP_array[q_WP_idx]   * q_WP_array[q_WP_idx])  /mu1(E_bound);
		double Eq_upper = 0.5*(q_WP_array[q_WP_idx+1] * q_WP_array[q_WP_idx+1])/mu1(E_bound);
		double E_com 	= 0.5*(Eq_upper + Eq_lower);
		double q 	 	= com_energy_to_com_q_momentum(E_com);
		double T_lab 	= com_momentum_to_lab_energy(q, E_bound);
		solve_config.T_lab_array[Tlab_idx] = T_lab;
	}

    solve_config.num_T_lab = num_T_lab;

	/* Calculate on-shell momentum and energy in centre-of-mass frame */
	solve_config.q_com_array = new double [num_T_lab];
	solve_config.E_com_array = new double [num_T_lab];
		
	for (size_t i=0; i<num_T_lab; i++){
		double q_com = lab_energy_to_com_momentum(solve_config.T_lab_array[i], E_bound);
		solve_config.q_com_array[i] = q_com;
		solve_config.E_com_array[i] = com_q_momentum_to_com_energy(q_com); 
	}

	printf("Locating on-shell q-momentum WP-indices for %zu on-shell energies ... \n", num_T_lab);
	solve_config.q_com_idx_array = new int [num_T_lab];
	for (size_t idx_Tlab=0; idx_Tlab<num_T_lab; idx_Tlab++){
		double q_com = solve_config.q_com_array[idx_Tlab];

		int idx_q_bin = -1;
		for (int q_idx_WP=0; q_idx_WP<Nq_WP; q_idx_WP++){
			if (q_WP_array[q_idx_WP]<q_com and q_com<q_WP_array[q_idx_WP+1]){
				idx_q_bin = q_idx_WP;
				break;
			}
		}

		if (idx_q_bin==-1){
			printf("On-shell kinetic energy Tlab=%.3f MeV doesn't exist in WP state space \n", solve_config.T_lab_array[idx_Tlab]);
			raise_error("Invalid Tlab entered. Exiting ...");
		}
		else{
			solve_config.q_com_idx_array[idx_Tlab] = idx_q_bin;
		}
	}
	printf(" - On-shell q-momentum WP bins found \n");
}

void find_deuteron_channels(solution_configuration& solve_config,
                            pw_3N_statespace pw_states){

	for (int chn_3N=0; chn_3N<pw_states.N_chn_3N; chn_3N++){
		/* Lower and upper limits on PW state space for channel */
		int idx_alpha_lower  = pw_states.chn_3N_idx_array[chn_3N];
		int idx_alpha_upper  = pw_states.chn_3N_idx_array[chn_3N+1];
		int Nalpha_in_3N_chn = idx_alpha_upper - idx_alpha_lower;

		/* Pointers to sub-arrays of PW state space corresponding to chn_3N */
		int* L_2N_subarray = &pw_states.L_2N_array[idx_alpha_lower];
		int* S_2N_subarray = &pw_states.S_2N_array[idx_alpha_lower];
		int* J_2N_subarray = &pw_states.J_2N_array[idx_alpha_lower];
		int* T_2N_subarray = &pw_states.T_2N_array[idx_alpha_lower];

		std::vector<int> deuteron_chn_indices;

		/* Find indices of deuteron-channels */
		for (int idx_alpha=0; idx_alpha<Nalpha_in_3N_chn; idx_alpha++){
			if (deuteron_L==L_2N_subarray[idx_alpha] &&
			    deuteron_S==S_2N_subarray[idx_alpha] &&
				deuteron_J==J_2N_subarray[idx_alpha] &&
				deuteron_T==T_2N_subarray[idx_alpha]){
				deuteron_chn_indices.push_back(idx_alpha);
			}
		}
			
		/* Copy vector content into array */
		int* chn_3N_deuteron_indices_array = new int [deuteron_chn_indices.size()];
		std::copy(deuteron_chn_indices.begin(), deuteron_chn_indices.end(), chn_3N_deuteron_indices_array);

		/* Add array length and pointer to book-keeping arrays */
		solve_config.deuteron_idx_arrays[chn_3N] = chn_3N_deuteron_indices_array;
		solve_config.deuteron_num_array [chn_3N] = deuteron_chn_indices.size();
	}
	printf(" - Nucleon-deuteron channel indices found \n");

	printf(" - Done \n");
}
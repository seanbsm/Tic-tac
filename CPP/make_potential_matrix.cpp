
#include "make_potential_matrix.h"

double extract_potential_element_from_array(int L, int Lp, int J, int S, bool coupled, double* V_array){
	if (J<abs(S-L) || J>S+L || J<abs(S-Lp) || J>S+Lp){
		std::cout << L <<" "<< Lp <<" "<< J <<" "<< S << std::endl;
		raise_error("Encountered unphysical state in LS-solver");
	}

	double potential_element = NAN;

	/* Sign-convention on off-diagonal coupled elements */
	double sgn = -1; // -1 USED IN BOUND STATE BENCHMARK

	if (coupled){
		if (L==Lp && L<J){         // --
			potential_element =  V_array[2];
		}
		else if (L==Lp && L>J){    // ++
			potential_element =  V_array[5];
		}
		else if (L<Lp){             // -+
			potential_element = sgn*V_array[4]; // 3 USED IN BOUND STATE BENCHMARK
		}
		else{                       // +-
			potential_element = sgn*V_array[3]; //4 USED IN BOUND STATE BENCHMARK
		}
	}
	else{
		if (L==J+1){         // 3P0 or (++)
			potential_element = V_array[5];
		}
		else if (L==J-1){			// (--)
			potential_element = V_array[2];
		}
		else if (S==0){             // S=0
			potential_element = V_array[0];
		}
		else if (S==1 && L==J){    // S=1
			potential_element = V_array[1];
		}
		else{
			potential_element = 0;
		}
	}

	return potential_element;
}

/* Construct 2N potential matrices <k|v|k_p> for all 3N partial wave channels */
void calculate_potential_matrices_array_in_WP_basis(double*  V_WP_unco_array, int num_2N_unco_states,
													double*  V_WP_coup_array, int num_2N_coup_states,
													fwp_statespace 	 fwp_states,
													pw_3N_statespace pw_states,
													potential_model* pot_ptr,
													run_params run_parameters){
	
	/* Make local pointers & variables for pw-statespace */
	int  Nalpha			= pw_states.Nalpha;
	int* L_2N_array		= pw_states.L_2N_array;
	int* S_2N_array		= pw_states.S_2N_array;
	int* J_2N_array		= pw_states.J_2N_array;
	int* T_2N_array		= pw_states.T_2N_array;
	int* two_T_3N_array	= pw_states.two_T_3N_array;
	/* Make local pointers & variables for FWP-statespace */
	int 	Np_WP		 = fwp_states.Np_WP;
	double* p_WP_array	 = fwp_states.p_WP_array;
	int 	Np_per_WP	 = fwp_states.Np_per_WP;
	double* p_array		 = fwp_states.p_array;
	double* wp_array	 = fwp_states.wp_array;
	double* norm_p_array = fwp_states.norm_p_array;
	double* fp_array 	 = fwp_states.fp_array;

	/* This test will be reused several times */
	bool tensor_force_true = (run_parameters.tensor_force==true);

	int J_2N_max = run_parameters.J_2N_max;

	int Tz_nn = -1;
	int Tz_np =  0;

	/* Potential-model input arrays */
	//double V_WP_elements [6];	// Isoscalar wave-packet potential elements (WP)
	//double V_IS_elements [6];	// Isoscalar (IS)
	//double V_nn_elements [6];	// neutron-neutron (nn)
	//double V_np_elements [6];	// neutron-proton (np)

	/* Potential matrix indexing */
	//int idx_V_WP_uncoupled   = 0;
	//int idx_V_WP_upper_left  = 0;
	//int idx_V_WP_upper_right = 0;
	//int idx_V_WP_lower_left  = 0;
	//int idx_V_WP_lower_right = 0;

	/* Temporary track-keeping arrays for avoiding repeat matrix-calculations (initialize to false) */
	int unco_array_size = num_2N_unco_states;
	int coup_array_size = num_2N_coup_states;
	
	bool* matrix_calculated_unco_array = new bool [unco_array_size];
	bool* matrix_calculated_coup_array = new bool [coup_array_size];
	for (int i=0; i<unco_array_size; i++){
		matrix_calculated_unco_array[i] = false;
	}
	for (int i=0; i<coup_array_size; i++){
		matrix_calculated_coup_array[i] = false;
	}

	printf("   - There are %d uncoupled 2N-channels and %d coupled 2N-channels \n", unco_array_size, coup_array_size);

	/* Row state */
	for (int idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
		int L_r = L_2N_array[idx_alpha_r];
		int S_r = S_2N_array[idx_alpha_r];
		int J_r = J_2N_array[idx_alpha_r];
		int T_r = T_2N_array[idx_alpha_r];

		int two_T_3N_r = two_T_3N_array[idx_alpha_r];
		
		/* Column state */
		for (int idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
			int L_c = L_2N_array[idx_alpha_c];
			int S_c = S_2N_array[idx_alpha_c];
			int J_c = J_2N_array[idx_alpha_c];
			int T_c = T_2N_array[idx_alpha_c];

			int two_T_3N_c = two_T_3N_array[idx_alpha_c];

			bool check_T = T_r==T_c;
			bool check_J = J_r==J_c;
			bool check_S = S_r==S_c;
			bool check_L = ( (tensor_force_true && abs(L_r-L_c)<=2) || L_r==L_c);

			/* Check if possible channel through interaction */
			if (check_T && check_J && check_S && check_L){

				/* Type of isoscalar interaction. The numbering is ordered in the same
				 * order as eq. (160) in Glockle et al., Phys. Rep. 274, 107 (1996) */
				int isoscalar_type = -1;

				/* Determine kind of isoscalar.
				 * For now the isospin-breaking terms have special treatment.
				 * Hopefully future code-iterations will make this completely general */
				bool state_1S0 = (S_r==0 && J_r==0 && L_r==0);
				if (state_1S0 && two_T_3N_r==3 && two_T_3N_c==3){
					isoscalar_type = 2; // Interaction can be either np or nn with IS symmetry-breaking
				}
				else if (state_1S0 && (two_T_3N_r==3 && two_T_3N_c==1) || (two_T_3N_r==1 && two_T_3N_c==3)){
					isoscalar_type = 3; // Interaction can be either np or nn with IS symmetry-conservation
				}
				else if (T_r==0 && two_T_3N_r==two_T_3N_c && two_T_3N_r==1){
					isoscalar_type = 0; // Interaction must be np
				}
				else if (T_r==1 && two_T_3N_r==two_T_3N_c && two_T_3N_r==1){
					isoscalar_type = 1; // Interaction can be either np or nn with IS symmetry-conservation
				}
				else{
					raise_error("Unknown isoscalar-potential encountered in potential-matrix construction!");
				}

				/* Detemine if this is a coupled channel.
				 * !!! With isospin symmetry-breaking we count 1S0 as a coupled matrix via T_3N-coupling !!! */
				bool coupled_matrix = false;
				bool coupled_model  = false;
				bool coupled_via_L_2N = (tensor_force_true && (L_r!=L_c || (L_r==L_c && L_r!=J_r && J_r!=0)));
				bool coupled_via_T_3N = (state_1S0==true && run_parameters.isospin_breaking_1S0==true);
				if (coupled_via_L_2N && coupled_via_T_3N){
					raise_error("Warning! Code has not been written to handle isospin-breaking in coupled channels!");
				}
				if (coupled_via_L_2N || coupled_via_T_3N){ // This counts 3P0 as uncoupled; used in matrix structure
					coupled_matrix  = true;
				}
				if (L_r!=L_c || (L_r==L_c && L_r!=J_r)){   // This counts 3P0 as coupled; used in potential models
					coupled_model   = true;
				}

				/* Unique 2N-channel indices */
				int chn_idx_V_coup = 0;
				int chn_idx_V_unco = 0;

				if (coupled_matrix){
					chn_idx_V_coup = unique_2N_idx(L_r, S_r, J_r, T_r, coupled_matrix, run_parameters);
				}
				else{
					chn_idx_V_unco = unique_2N_idx(L_r, S_r, J_r, T_r, coupled_matrix, run_parameters);
				}

				/* Check if calculation has already been performed in some other alpha'-alpha interaction */
				if (coupled_matrix){
					if (matrix_calculated_coup_array[chn_idx_V_coup]==true){
						continue;
					}
					else{
						matrix_calculated_coup_array[chn_idx_V_coup] = true;
					}
				}
				else{
					if (matrix_calculated_unco_array[chn_idx_V_unco]==true){
						continue;
					}
					else{
						matrix_calculated_unco_array[chn_idx_V_unco] = true;
					}
				}

				printf("     - Working on matrix <L'=%d|V(S=%d, J=%d, T=%d)|L=%d> ", L_r, S_r, J_r, T_r, L_c);
				if (coupled_via_L_2N){
					printf("(coupled matrix via L_2N -> calculation includes all couplings of L' and L.)");
				}
				if (coupled_via_T_3N){
					printf("(coupled matrix via T_3N -> isospin-breaking.)");
				}
				printf("\n");

				/* Pre-calculate matrices for efficient parameter sampling (ONLY USED IF FIRST TIME CALLED)
				 * for current NN-channel */
				//int Np = 0;
				//if (run_parameters.midpoint_approx==true){
				//	Np = Np_WP;
				//}
				//else{
				//	Np = Np_WP*Np_per_WP;
				//}
				//pot_ptr->setup_store_matrices(p_array, Np, coupled_model, S_r, J_r, T_r, Tz_np);
				//if (isoscalar_type>=1){
				//	pot_ptr->setup_store_matrices(p_array, Np, coupled_model, S_r, J_r, T_r, Tz_nn);
				//}

				/* Row p-momentum index loop */
				//#pragma omp parallel	// Parallel looping is not thread-safe for chiral potential or funtional with nijmegen potential - out of my control...
				//{
				/* Potential-model input arrays */
				double V_WP_elements [6];
				double V_IS_elements [6];
				double V_nn_elements [6];
				double V_np_elements [6];
				/* Potential matrix indexing */
				int idx_V_WP_uncoupled   = 0;
				int idx_V_WP_upper_left  = 0;
				int idx_V_WP_upper_right = 0;
				int idx_V_WP_lower_left  = 0;
				int idx_V_WP_lower_right = 0;
				//#pragma for
				for (int idx_bin_r=0; idx_bin_r<Np_WP; idx_bin_r++){
					/* WP normalization */
					double N_r = norm_p_array[idx_bin_r];

					/* Column p-momentum index loop */
					for (int idx_bin_c=0; idx_bin_c<Np_WP; idx_bin_c++){
						/* WP normalization */
						double N_c = norm_p_array[idx_bin_c];

						/* Potential matrix indexing
						 * Indexing format: (channel index)*(num rows)*(num columns) + (row index)*(row length) + (column index) */
						if (coupled_matrix){
							idx_V_WP_upper_left  = chn_idx_V_coup*4*Np_WP*Np_WP +  idx_bin_r         *2*Np_WP + idx_bin_c;
							idx_V_WP_upper_right = chn_idx_V_coup*4*Np_WP*Np_WP +  idx_bin_r         *2*Np_WP + idx_bin_c + Np_WP;
							idx_V_WP_lower_left  = chn_idx_V_coup*4*Np_WP*Np_WP + (idx_bin_r + Np_WP)*2*Np_WP + idx_bin_c;
							idx_V_WP_lower_right = chn_idx_V_coup*4*Np_WP*Np_WP + (idx_bin_r + Np_WP)*2*Np_WP + idx_bin_c + Np_WP;
						
						}
						else{
							idx_V_WP_uncoupled   = chn_idx_V_unco*Np_WP*Np_WP + idx_bin_r*Np_WP + idx_bin_c;
						}
						
						/* Reset quadrature summation array V_WP_elements */
						for (int idx_element=0; idx_element<6; idx_element++){
							V_WP_elements[idx_element] = 0;
						}

						/* Calculate matrix element, using either mid-point approximation or quadrature */
						if (run_parameters.midpoint_approx==true){
							
							/* Average momenta */
							double p_in  = p_array[idx_bin_c];
							double p_out = p_array[idx_bin_r];

							/* Weigth-function terms */
							double f_in  = fp_array[idx_bin_c];
							double f_out = fp_array[idx_bin_r];

							/* Momentum bin-width */
							double d_c   = p_WP_array[idx_bin_c + 1] - p_WP_array[idx_bin_c];
							double d_r   = p_WP_array[idx_bin_r + 1] - p_WP_array[idx_bin_r];

							/* Row and column indices of p_out and p_in, respectively */
							int i = idx_bin_r;
							int j = idx_bin_c;

							/* We create an isoscalar potential */
							if (isoscalar_type==0){ // Interaction must be np
								pot_ptr->V(i, j, p_in, p_out, coupled_model, S_r, J_r, T_r, Tz_np, V_IS_elements);
							}
							else if (isoscalar_type==1){ // Interaction can be either nn or np with IS symmetry-conservation
								pot_ptr->V(i, j, p_in, p_out, coupled_model, S_r, J_r, T_r, Tz_nn, V_nn_elements);
								pot_ptr->V(i, j, p_in, p_out, coupled_model, S_r, J_r, T_r, Tz_np, V_np_elements);
								for (int idx_element=0; idx_element<6; idx_element++){
									V_IS_elements[idx_element] = (2./3)*V_nn_elements[idx_element] + (1./3)*V_np_elements[idx_element];
								}
							}
							else if (isoscalar_type==2 || isoscalar_type==3){ // Interaction can be either np or nn with IS symmetry-breaking
								pot_ptr->V(i, j, p_in, p_out, coupled_model, S_r, J_r, T_r, Tz_nn, V_nn_elements);
								pot_ptr->V(i, j, p_in, p_out, coupled_model, S_r, J_r, T_r, Tz_np, V_np_elements);
								
								/* Extract 1S0-element, using coupling as according to tensor-force only (i.e. no coupling for 1S0) */
								double V_nn_1S0_element = extract_potential_element_from_array(L_r, L_c, J_r, S_r, false, V_nn_elements);
								double V_np_1S0_element = extract_potential_element_from_array(L_r, L_c, J_r, S_r, false, V_np_elements);

								/* Define our own, convenient ordering of V_IS_elements:
								 * idx 2: top left block		(isoscalar type 1)
								 * idx 3: top right block		(isoscalar type 2)
								 * idx 4: bottom left block		(isoscalar type 2)
								 * idx 5: bottom right block	(isoscalar type 3) */
								V_IS_elements[2] = 		 (2./3)*  V_nn_1S0_element + (1./3)*V_np_1S0_element;
								V_IS_elements[3] = (sqrt(2)/3.)* (V_nn_1S0_element - 	    V_np_1S0_element);
								V_IS_elements[4] = (sqrt(2)/3.)* (V_nn_1S0_element - 	    V_np_1S0_element);
								V_IS_elements[5] = 		 (1./3)*  V_nn_1S0_element + (2./3)*V_np_1S0_element;
								
							}
							else{
								raise_error("Unknown isoscalar potential encountered.");
							}
							
							/* We integrate into wave-packet potential and normalize */
							for (int idx_element=0; idx_element<6; idx_element++){
								V_WP_elements[idx_element] += p_in*p_out *f_in*f_out *d_r*d_c * V_IS_elements[idx_element] / (N_r*N_c);
							}
						}
						else{
							/* Loop over quadrature inside row and column cells */
							for (int idx_p_r=0; idx_p_r<Np_per_WP; idx_p_r++){

								/* Integration setup */
								double  p_r     =  p_array[idx_bin_r*Np_per_WP + idx_p_r];
								double wp_r     = wp_array[idx_bin_r*Np_per_WP + idx_p_r];
								double wp_p_f_r = fp_array[idx_bin_r*Np_per_WP + idx_p_r] * wp_r*p_r;
								
								/* Variable change for easier readablity: <p_out|v|p_in> */
								double p_out = p_r;

								/* Row index of p_out */
								int i = idx_bin_r*Np_per_WP + idx_p_r;

								for (int idx_p_c=0; idx_p_c<Np_per_WP; idx_p_c++){

									/* Integration setup */
									double  p_c 	=  p_array[idx_bin_c*Np_per_WP + idx_p_c];
									double wp_c 	= wp_array[idx_bin_c*Np_per_WP + idx_p_c];
									double wp_p_f_c = fp_array[idx_bin_c*Np_per_WP + idx_p_c] * wp_c*p_c;

									/* Variable change for easier readablity: <p_out|v|p_in> */
									double p_in = p_c;

									/* Column index of p_in */
									int j = idx_bin_c*Np_per_WP + idx_p_c;

									/* Prefactors for quadrature and WP-weighting */
									double integral_factors = wp_p_f_r * wp_p_f_c;

									/* We create an isoscalar potential */
									if (coupled_via_T_3N==false && isoscalar_type==0){ // Interaction must be np
										pot_ptr->V(i, j, p_in, p_out, coupled_model, S_r, J_r, T_r, Tz_np, V_IS_elements);
									}
									else if (coupled_via_T_3N==false && isoscalar_type==1){ // Interaction can be either nn or np with IS symmetry-conservation
										pot_ptr->V(i, j, p_in, p_out, coupled_model, S_r, J_r, T_r, Tz_nn, V_nn_elements);
										pot_ptr->V(i, j, p_in, p_out, coupled_model, S_r, J_r, T_r, Tz_np, V_np_elements);
										for (int idx_element=0; idx_element<6; idx_element++){
											V_IS_elements[idx_element] = (2./3)*V_nn_elements[idx_element] + (1./3)*V_np_elements[idx_element];
										}
									}
									else if (coupled_via_T_3N==true){ // Interaction can be either np or nn with IS symmetry-breaking
										pot_ptr->V(i, j, p_in, p_out, coupled_model, S_r, J_r, T_r, Tz_nn, V_nn_elements);
										pot_ptr->V(i, j, p_in, p_out, coupled_model, S_r, J_r, T_r, Tz_np, V_np_elements);

										/* Extract 1S0-element, using coupling as according to tensor-force only (i.e. no coupling for 1S0) */
										double V_nn_1S0_element = extract_potential_element_from_array(L_r, L_c, J_r, S_r, false, V_nn_elements);
										double V_np_1S0_element = extract_potential_element_from_array(L_r, L_c, J_r, S_r, false, V_np_elements);

										/* Define our own, convenient ordering of V_IS_elements:
										 * idx 2: top left block		(isoscalar type 1)
										 * idx 3: top right block		(isoscalar type 2)
										 * idx 4: bottom left block		(isoscalar type 2)
										 * idx 5: bottom right block	(isoscalar type 3) */
										V_IS_elements[2] = 		 (2./3)*  V_nn_1S0_element + (1./3)*V_np_1S0_element;
										V_IS_elements[3] = (sqrt(2)/3.)* (V_nn_1S0_element - 	    V_np_1S0_element);
										V_IS_elements[4] = (sqrt(2)/3.)* (V_nn_1S0_element - 	    V_np_1S0_element);
										V_IS_elements[5] = 		 (1./3)*  V_nn_1S0_element + (2./3)*V_np_1S0_element;
									}
									else{
										raise_error("Unknown isoscalar potential encountered.");
									}

									/* We integrate into wave-packet potential and normalize */
									for (int idx_element=0; idx_element<6; idx_element++){
										V_WP_elements[idx_element] += integral_factors * V_IS_elements[idx_element]/(N_r*N_c);
									}
								}
							}
						}
						
						/* Write element to potential matrix V_array */
						if (coupled_matrix){
							if (coupled_via_L_2N){
								V_WP_coup_array[idx_V_WP_upper_left]  = extract_potential_element_from_array(J_r-1, J_r-1, J_r, S_r, coupled_matrix, V_WP_elements);
								V_WP_coup_array[idx_V_WP_upper_right] = extract_potential_element_from_array(J_r-1, J_r+1, J_r, S_r, coupled_matrix, V_WP_elements);
								V_WP_coup_array[idx_V_WP_lower_left]  = extract_potential_element_from_array(J_r+1, J_r-1, J_r, S_r, coupled_matrix, V_WP_elements);
								V_WP_coup_array[idx_V_WP_lower_right] = extract_potential_element_from_array(J_r+1, J_r+1, J_r, S_r, coupled_matrix, V_WP_elements);
							}
							else if (coupled_via_T_3N){
								V_WP_coup_array[idx_V_WP_upper_left]  = V_WP_elements[2];
								V_WP_coup_array[idx_V_WP_upper_right] = V_WP_elements[3];
								V_WP_coup_array[idx_V_WP_lower_left]  = V_WP_elements[4];
								V_WP_coup_array[idx_V_WP_lower_right] = V_WP_elements[5];
							}
							//if (J_r==1 and L_r==0){
							//    std::cout << idx_bin_r <<" "<< idx_bin_c <<" "<< V_WP_coup_array[idx_V_WP_upper_left] << std::endl;
							//    std::cout << idx_bin_r <<" "<< idx_bin_c <<" "<< V_WP_coup_array[idx_V_WP_upper_right] << std::endl;
							//    std::cout << idx_bin_r <<" "<< idx_bin_c <<" "<< V_WP_coup_array[idx_V_WP_lower_left] << std::endl;
							//    std::cout << idx_bin_r <<" "<< idx_bin_c <<" "<< V_WP_coup_array[idx_V_WP_lower_right] << std::endl;
							//}
						}
						else{
							V_WP_unco_array[idx_V_WP_uncoupled]   = extract_potential_element_from_array(L_r, L_c, J_r, S_r, coupled_matrix, V_WP_elements);
						}
					}
				}
				//}
			}
		}
	}
}
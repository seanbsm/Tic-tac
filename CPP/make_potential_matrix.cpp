
#include "make_potential_matrix.h"

double extract_potential_element_from_array(int L, int Lp, int J, int S, bool coupled, double* V_array){
	if (J<abs(S-L) or J>S+L or J<abs(S-Lp) or J>S+Lp){
		std::cout << L <<" "<< Lp <<" "<< J <<" "<< S << std::endl;
		raise_error("Encountered unphysical state in LS-solver");
	}

	double potential_element = NAN;

	/* Sign-convention on off-diagonal coupled elements */
	double sgn = -1; // -1 USED IN BOUND STATE BENCHMARK

	if (coupled){
		if (L==Lp and L<J){         // --
			potential_element =  V_array[2];
		}
		else if (L==Lp and L>J){    // ++
			potential_element =  V_array[5];
		}
		else if (L<Lp){             // -+
			potential_element = sgn*V_array[3]; // 3 USED IN BOUND STATE BENCHMARK
		}
		else{                       // +-
			potential_element = sgn*V_array[4]; //4 USED IN BOUND STATE BENCHMARK
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
		else if (S==1 and L==J){    // S=1
			potential_element = V_array[1];
		}
		else{
			potential_element = 0;
		}
	}

	return potential_element;
}

/* Construct 2N potential matrices <k|v|k_p> for all 3N partial wave channels */
void calculate_potential_matrices_array_in_WP_basis(double*  V_WP_unco_array,
													double*  V_WP_coup_array,
													bool tensor_force_true,
													bool mid_point_approximation,
													int Np_WP, double* p_WP_array,
													int Np_per_WP, double* p_array, double* wp_array,
													int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array,
													int J_2N_max,
													potential_model* pot_ptr_nn,
													potential_model* pot_ptr_np){
	
	/* Potential-model input arrays */
	double V_WP_elements [6];	// Isoscalar wave-packet potential elements (WP)
	double V_IS_elements [6];	// Isoscalar (IS)
	double V_nn_elements [6];	// neutron-neutron (nn)
	double V_np_elements [6];	// neutron-proton (np)

	/* Potential matrix indexing */
	int idx_V_WP_uncoupled   = 0;
	int idx_V_WP_upper_left  = 0;
	int idx_V_WP_upper_right = 0;
	int idx_V_WP_lower_left  = 0;
	int idx_V_WP_lower_right = 0;

	/* Temporary track-keeping arrays for avoiding repeat matrix-calculations (initialize to false) */
	int unco_array_size = 0;
	int coup_array_size = 0;
	if (tensor_force_true){
		unco_array_size = 2*(J_2N_max+1);
		coup_array_size =    J_2N_max;
	}
	else{
		unco_array_size = 4*J_2N_max + 2;
	}
	bool* matrix_calculated_unco_array = new bool [unco_array_size];
	bool* matrix_calculated_coup_array = new bool [coup_array_size];
	for (int i=0; i<unco_array_size; i++){
		matrix_calculated_unco_array[i] = false;
	}
	for (int i=0; i<coup_array_size; i++){
		matrix_calculated_coup_array[i] = false;
	}

	printf(" - There are %d uncoupled 2N-channels and %d coupled 2N-channels \n", unco_array_size, coup_array_size);

	/* Row state */
	for (int idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
		int L_r = L_2N_array[idx_alpha_r];
		int S_r = S_2N_array[idx_alpha_r];
		int J_r = J_2N_array[idx_alpha_r];
		int T_r = T_2N_array[idx_alpha_r];
		
		/* Column state */
		for (int idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
			int L_c = L_2N_array[idx_alpha_c];
			int S_c = S_2N_array[idx_alpha_c];
			int J_c = J_2N_array[idx_alpha_c];
			int T_c = T_2N_array[idx_alpha_c];
			

			bool check_T = T_r==T_c;
			bool check_J = J_r==J_c;
			bool check_S = S_r==S_c;
			bool check_L = ( (tensor_force_true && abs(L_r-L_c)<=2) || L_r==L_c);


			/* Check if possible channel through interaction */
			if (check_T && check_J && check_S && check_L){

				/* Detemine if this is a coupled channel */
				bool coupled_matrix = false;
				bool coupled_model  = false;
				if (tensor_force_true && (L_r!=L_c or (L_r==L_c and L_r!=J_r and J_r!=0)) ){    // This counts 3P0 as uncoupled; used in matrix structure
					coupled_matrix  = true;
				}
				if (L_r!=L_c or (L_r==L_c and L_r!=J_r)){               // This counts 3P0 as coupled; used in potential models
					coupled_model   = true;
				}

				/* Unique 2N-channel indices */
				int chn_idx_V_coup = 0;
				int chn_idx_V_unco = 0;
				if (coupled_matrix){
					chn_idx_V_coup = unique_2N_idx(L_r, S_r, J_r, T_r, tensor_force_true, coupled_matrix);
				}
				else{
					chn_idx_V_unco = unique_2N_idx(L_r, S_r, J_r, T_r, tensor_force_true, coupled_matrix);
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

				printf("   - Working on matrix <L'=%d|V(S=%d, J=%d, T=%d)|L=%d> ", L_r, S_r, J_r, T_r, L_c);
				if (coupled_matrix){
					printf("(coupled matrix -> calculation includes all couplings of L' and L.)");
				}
				printf("\n");

				/* Skip redundant calculations by only doing the coupled calculation when L_r<L_c */
				//if (coupled_matrix){
				//	if ( (L_r<L_c)==false ){
				//		continue;
				//	}
				//}

				/* Row p-momentum index loop */
				for (int idx_bin_r=0; idx_bin_r<Np_WP; idx_bin_r++){
					double bin_r_lower = p_WP_array[idx_bin_r];
					double bin_r_upper = p_WP_array[idx_bin_r + 1];
					double N_r = p_normalization(bin_r_lower, bin_r_upper);

					double*  p_array_ptr_r =  &p_array[idx_bin_r*Np_per_WP];
					double* wp_array_ptr_r = &wp_array[idx_bin_r*Np_per_WP];

					/* Column p-momentum index loop */
					for (int idx_bin_c=0; idx_bin_c<Np_WP; idx_bin_c++){
						double bin_c_lower = p_WP_array[idx_bin_c];
						double bin_c_upper = p_WP_array[idx_bin_c + 1];
						double N_c = p_normalization(bin_c_lower, bin_c_upper);

						double*  p_array_ptr_c =  &p_array[idx_bin_c*Np_per_WP];
						double* wp_array_ptr_c = &wp_array[idx_bin_c*Np_per_WP];

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
						if (mid_point_approximation){
							/* Average momentum */
							double p_in  = 0.5*(bin_c_lower + bin_c_upper);
							double p_out = 0.5*(bin_r_lower + bin_r_upper);

							/* Weigth-function terms */
							double f_in  = p_weight_function(p_in);
							double f_out = p_weight_function(p_out);

							/* Momentum bin-width */
							double d_c   = bin_c_upper - bin_c_lower;
							double d_r   = bin_r_upper - bin_r_lower;

							/* We create an isoscalar potential */
							if (T_r==1){ // Interaction can be either nn or np
								pot_ptr_nn->V(p_in, p_out, coupled_model, S_r, J_r, T_r, V_nn_elements);
								pot_ptr_np->V(p_in, p_out, coupled_model, S_r, J_r, T_r, V_np_elements);
								for (int idx_element=0; idx_element<6; idx_element++){
									V_IS_elements[idx_element] = (1./3)*V_np_elements[idx_element] + (2./3)*V_nn_elements[idx_element];
								}
							}
							else{ 	   // Interaction must be np
								pot_ptr_np->V(p_in, p_out, coupled_model, S_r, J_r, T_r, V_IS_elements);
							}
							
							/* We integrate into wave-packet potential and normalize */
							for (int idx_element=0; idx_element<6; idx_element++){
								//V_WP_elements[idx_element] += p_in*p_out *f_in*f_out *d_r*d_c * V_IS_elements[idx_element] / (N_r*N_c);

								double p_r_2 = bin_r_upper*bin_r_upper - bin_r_lower*bin_r_lower;
								double p_c_2 = bin_c_upper*bin_c_upper - bin_c_lower*bin_c_lower;
								V_WP_elements[idx_element] += p_r_2 * p_c_2 * V_IS_elements[idx_element] / (4*N_r*N_c);
							}
						}
						else{
							/* Loop over quadrature inside row and column cells */
							for (int idx_p_r=0; idx_p_r<Np_per_WP; idx_p_r++){
								double  p_r  =  p_array_ptr_r[idx_p_r];
								double wp_r  = wp_array_ptr_r[idx_p_r];
								double p_out = p_r;     // variable change for easier readablity: <p_out|v|p_in>

								for (int idx_p_c=0; idx_p_c<Np_per_WP; idx_p_c++){
									double  p_c =  p_array_ptr_c[idx_p_c];
									double wp_c = wp_array_ptr_c[idx_p_c];
									double p_in = p_c;  // variable change for easier readablity: <p_out|v|p_in>

									double wp_p_f_r = wp_r*p_r*p_weight_function(p_r);
									double wp_p_f_c = wp_c*p_c*p_weight_function(p_c);
									double integral_factors = wp_p_f_r * wp_p_f_c;

									/* We create an isoscalar potential */
									if (T_r==1){ // Interaction can be either nn or np
										pot_ptr_nn->V(p_in, p_out, coupled_model, S_r, J_r, T_r, V_nn_elements);
										pot_ptr_np->V(p_in, p_out, coupled_model, S_r, J_r, T_r, V_np_elements);
										for (int idx_element=0; idx_element<6; idx_element++){
											V_IS_elements[idx_element] = (1./3)*V_np_elements[idx_element] + (2./3)*V_nn_elements[idx_element];
										}
									}
									else{ 	   // Interaction must be np
										pot_ptr_np->V(p_in, p_out, coupled_model, S_r, J_r, T_r, V_IS_elements);
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
							V_WP_coup_array[idx_V_WP_upper_left]  = extract_potential_element_from_array(J_r-1, J_r-1, J_r, S_r, coupled_matrix, V_WP_elements);
							V_WP_coup_array[idx_V_WP_upper_right] = extract_potential_element_from_array(J_r-1, J_r+1, J_r, S_r, coupled_matrix, V_WP_elements);
							V_WP_coup_array[idx_V_WP_lower_left]  = extract_potential_element_from_array(J_r+1, J_r-1, J_r, S_r, coupled_matrix, V_WP_elements);
							V_WP_coup_array[idx_V_WP_lower_right] = extract_potential_element_from_array(J_r+1, J_r+1, J_r, S_r, coupled_matrix, V_WP_elements);
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
			}
		}
	}
}
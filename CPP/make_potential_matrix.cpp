
#include "make_potential_matrix.h"

double extract_potential_element_from_array(int L, int Lp, int J, int S, bool coupled, double* V_array){
	if (J<abs(S-L) or J>S+L or J<abs(S-Lp) or J>S+Lp){
		std::cout << L <<" "<< Lp <<" "<< J <<" "<< S << std::endl;
		raise_error("Encountered unphysical state in LS-solver");
	}

    double potential_element = NAN;

    if (coupled){
        if (L==Lp and L<J){         // --
            potential_element = V_array[2];
        }
        else if (L==Lp and L>J){    // ++
            potential_element = V_array[5];
        }
        else if (L<Lp){             // +-
            potential_element = V_array[3];
        }
        else{                       // -+
            potential_element = V_array[4];
        }
    }
    else{
        if (J==0 and L!=J){         // 3P0 (++)
            potential_element = V_array[5];
        }
        else if (S==0){             // S=0
            potential_element = V_array[0];
        }
        else{                       // S=1
            potential_element = V_array[1];
        }
    }

    return potential_element;
}

/* Construct 2N potential matrices <k|v|k_p> for all 3N partial wave channels */
void calculate_potential_matrices_array(double* V_WP_unco_array,
                                        double* V_WP_coup_array,
                                        int Np_WP, double* p_WP_array,
                                        int Np, double* p_array, double* wp_array,
                                        int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array,
                                        potential_model* pot_ptr_nn,
                                        potential_model* pot_ptr_np){
    
    double V_WP_elements [6];	// Isoscalar wave-packet potential elements (WP)
    double V_IS_elements [6];	// Isoscalar (IS)
	double V_nn_elements [6];	// neutron-neutron (nn)
	double V_np_elements [6];	// neutron-proton (np)

    /* Potential matrix indexing */
    int idx_V_WP_uncoupled   = 0;
    int idx_V_WP_upper_left  = 0;
	int idx_V_WP_Upper_right = 0;
	int idx_V_WP_lower_left  = 0;
	int idx_V_WP_lower_right = 0;

	double p_r=0, p_in=0, p_c=0, p_out=0;
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

            /* Check if possible channel through interaction */
            if (T_r==T_c and J_r==J_c and S_r==S_c and abs(L_r-L_c)<=2){
                /* Detemine if this is a coupled channel */
	            bool coupled = false;
	            if (L_r!=L_c or (L_r==L_c and L_r!=J_r and J_r!=0)){
	            	coupled  = true;
	            }

                /* Skip redundant calculations by only doing the coupled calculation when L_r<L_c */
                if (coupled){
                    if ( (L_r<L_c)==false ){
                        continue;
                    }
                }

                /* Row index loop */
                for (int idx_bin_r=0; idx_bin_r<Np_WP; idx_bin_r++){
                    double bin_r_lower = p_WP_array[idx_bin_r];
                    double bin_r_upper = p_WP_array[idx_bin_r + 1];
                    double N_r = p_normalisation(bin_r_lower, bin_r_upper);
                    
                    /* Column index loop */
                    for (int idx_bin_c=0; idx_bin_c<Np_WP; idx_bin_c++){
                        double bin_c_lower = p_WP_array[idx_bin_c];
                        double bin_c_upper = p_WP_array[idx_bin_c + 1];
                        double N_c = p_normalisation(bin_c_lower, bin_c_upper);

                        /* Potential matrix indexing */
                        if (coupled){
                            int step_V_coup      = J_r-1;
			            	idx_V_WP_upper_left  = step_V_coup*4*Np_WP*Np_WP +  idx_bin_r         *2*Np_WP + idx_bin_c;
			            	idx_V_WP_Upper_right = step_V_coup*4*Np_WP*Np_WP +  idx_bin_r         *2*Np_WP + idx_bin_c + Np_WP;
			            	idx_V_WP_lower_left  = step_V_coup*4*Np_WP*Np_WP + (idx_bin_r + Np_WP)*2*Np_WP + idx_bin_c;
			            	idx_V_WP_lower_right = step_V_coup*4*Np_WP*Np_WP + (idx_bin_r + Np_WP)*2*Np_WP + idx_bin_c + Np_WP;
                        
                        }
			            else{
                            int step_V_unco      = L_r + S_r + (J_r!=0) - (J_r==0 and L_r!=J_r);   // This indexing gives room for the 3P0-wave
			            	idx_V_WP_uncoupled   = step_V_unco*Np_WP*Np_WP + idx_bin_r*Np_WP + idx_bin_c;
                        }

                        /* Reset quadrature summation array V_WP_elements */
	                    for (int idx_element=0; idx_element<6; idx_element++){
	                    	V_WP_elements[idx_element] = 0;
                        }

                        /* Loop over quadrature inside row and column bins */
                        for (int idx_p_r=0; idx_p_r<Np; idx_p_r++){
                            double p_r   = p_array[idx_p_r];
                            double wp_r  = wp_array[idx_p_r];
                            double p_out = p_r; //*hbarc;

                            for (int idx_p_c=0; idx_p_c<Np; idx_p_c++){
                                double p_c  = p_array[idx_p_c];
                                double wp_c = wp_array[idx_p_c];
                                double p_in = p_c; //*hbarc;

	                            /* We create an isoscalar potential */
	                            if (T_r==1){ // Interaction can be either nn or np
                                    pot_ptr_nn->V(p_in, p_out, coupled, S_r, J_r, T_r, V_nn_elements);
                                    pot_ptr_np->V(p_in, p_out, coupled, S_r, J_r, T_r, V_np_elements);
	                            	for (int idx_element=0; idx_element<6; idx_element++){
	                            		V_IS_elements[idx_element] = (1./3)*V_np_elements[idx_element] + (2./3)*V_nn_elements[idx_element];
	                            	}
                                }
                                else{ 	   // Interaction must be np
                                    pot_ptr_np->V(p_in, p_out, coupled, S_r, J_r, T_r, V_IS_elements);
                                }

                                /* We integrate into wave-packet potential */
	                            for (int idx_element=0; idx_element<6; idx_element++){
	                            	V_WP_elements[idx_element] += p_r*p_c * wp_r*wp_c * V_IS_elements[idx_element]/(sqrt(N_r*N_c));
                                }
                            }
                        }

	                    /* Write element to potential matrix V_array */
                        if (coupled){
                            V_WP_coup_array[idx_V_WP_upper_left]  = extract_potential_element_from_array(J_r-1, J_r-1, J_r, S_r, coupled, V_WP_elements);
			            	V_WP_coup_array[idx_V_WP_Upper_right] = extract_potential_element_from_array(J_r-1, J_r+1, J_r, S_r, coupled, V_WP_elements);
			            	V_WP_coup_array[idx_V_WP_lower_left]  = extract_potential_element_from_array(J_r+1, J_r-1, J_r, S_r, coupled, V_WP_elements);
			            	V_WP_coup_array[idx_V_WP_lower_right] = extract_potential_element_from_array(J_r+1, J_r+1, J_r, S_r, coupled, V_WP_elements);
                        
                        }
			            else{
                            V_WP_unco_array[idx_V_WP_uncoupled] = extract_potential_element_from_array(L_r, L_c, J_r, S_r, coupled, V_WP_elements);

                        }
                    }
                }
            }
        }
    }
}

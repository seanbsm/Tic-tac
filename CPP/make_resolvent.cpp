
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
cdouble resolvent_bound_continuum(double E,
                                  double q_bin_upper,
                                  double q_bin_lower,
                                  double Eb){
    
    /* Bin boundaries in energy */
    double Eq_lower = q_bin_lower * q_bin_lower / MN;
    double Eq_upper = q_bin_upper * q_bin_upper / MN;

    /* Temporary variable */
    double e = Eb - E;

    /* Energy width of SWPs (D is short for Delta) for q momenta */
    double Dq = Eq_upper - Eq_lower;

    /* Real part of the BC resolvent */
    double Re_R = log( abs( (Eq_lower + e)/(Eq_upper + e) ) ) / (Dq); 
    
    /* Imaginary part of the BC resolvent */
    double Im_R = (  heaviside_step_function( Eq_upper + e )
                   - heaviside_step_function( Eq_lower + e ) ) * (-M_PI)/ (Dq); 
    
    /* Return the complex BC resolvent term Q */
    return {Re_R, Im_R};
}

/* See header-file, commentary (A), for explanation of notation and equations */
cdouble resolvent_continuum_continuum(double E,
                                      double q_bin_upper,
                                      double q_bin_lower,
                                      double p_bin_upper,
                                      double p_bin_lower){
    
    /* Bin boundaries in energy */
    double Eq_lower = q_bin_lower * q_bin_lower / MN;
    double Eq_upper = q_bin_upper * q_bin_upper / MN;

    double Ep_lower = p_bin_lower * p_bin_lower / MN;
    double Ep_upper = p_bin_upper * p_bin_upper / MN;

    /* Energy widths of SWPs (D is short for Delta) for p and q momenta */
    double Dq = Eq_upper - Eq_lower;
    double Dp = Ep_upper - Ep_lower;

    /* Temporary variables */
    double DM = 0.5*(Dp - Dq);
    double DP = 0.5*(Dp + Dq);
    double D  = Ep_upper + Eq_upper - E;

    /* Real part of the CC resolvent */
    double Re_Q = (  (D+DM) * log( abs(D+DM) )
                   + (D-DM) * log( abs(D-DM) )
                   - (D+DP) * log( abs(D+DP) )
                   - (D-DP) * log( abs(D-DP) ) ) / (Dp*Dq); 
    
    /* Imaginary part of the CC resolvent */
    double Im_Q = (  (D+DM) * heaviside_step_function( D+DM )
                   + (D-DM) * heaviside_step_function( D-DM )
                   - (D+DP) * heaviside_step_function( D+DP )
                   - (D-DP) * heaviside_step_function( D-DP ) ) * (-M_PI)/ (Dp*Dq); 
    
    /* Return the complex CC resolvent term Q */
    return {Re_Q, Im_Q};
}

void calculate_resolvent_array_in_SWP_basis(cdouble** G_array,
                                            double  E,
                                            int     Np_WP, double* p_SWP_unco_array, double* p_SWP_coup_array,
					                        int     Nq_WP, double* q_WP_array,
                                            int     N_chn_3N, int* chn_3N_idx_array,
					                        int     Nalpha,
                                            int*    L_2N_array,
                                            int*    S_2N_array,
                                            int*    J_2N_array,
                                            int*    T_2N_array){

    /* Pointer to either p_SWP_unco_array or p_SWP_coup_array,
     * which is determined by whether the channel is coupled or not */
    double* p_SWP_array_ptr = NULL;

    for (int chn_3N=0; chn_3N<N_chn_3N; chn_3N++){
        /* Partial-waves in current 3N channel */
        int alpha_idx_lower = chn_3N_idx_array[chn_3N];
        int alpha_idx_upper = chn_3N_idx_array[chn_3N];
        
        /* Note that chn_3N_idx_array[N_chn_3N] = Nalpha */
        int Nalpha_block = alpha_idx_upper - alpha_idx_lower;

        /* Allocate G_array sub-array */
        G_array[chn_3N] = new cdouble [Nalpha_block * Nq_WP * Np_WP];

        /* Create pointer to sub-array for simplicity */
        cdouble* G_subarray = G_array[chn_3N];

        /* Loop over states along resolvent diagonal */
        for (int idx_alpha=alpha_idx_lower; idx_alpha<alpha_idx_upper; idx_alpha++){
    
            int L = L_2N_array[idx_alpha];
            int S = S_2N_array[idx_alpha];
            int J = J_2N_array[idx_alpha];
            int T = T_2N_array[idx_alpha];

            /* Detemine if this is a coupled channel */
	        if (L!=J and J!=0){ // This counts 3P0 as uncoupled (which is intended)
	        	int mat_dim = 2*Np_WP;
                int chn_2N_idx = J-1;
                p_SWP_array_ptr = &p_SWP_coup_array[chn_2N_idx * mat_dim];
	        }
            else{
                int mat_dim = Np_WP;
                int chn_2N_idx = 2*J + S;
                p_SWP_array_ptr = &p_SWP_unco_array[chn_2N_idx * mat_dim];
            }

            /* p-momentum index loop */
            for (int idx_p_bin=0; idx_p_bin<Np_WP; idx_p_bin++){

                /* Upper and lower boundaries of current p-bin */
                double p_bin_lower = p_SWP_array_ptr[2*idx_p_bin    ];
                double p_bin_upper = p_SWP_array_ptr[2*idx_p_bin + 1];

                /* Bound state check, given by p_bin_lower if bound state exists */
                bool bound_state_exists = false;
                double Eb = 0;
                if (p_bin_lower<0){
                    Eb = p_bin_lower;
                    bound_state_exists = true;
                }

                /* q-momentum index loop */            
                for (int idx_q_bin=0; idx_q_bin<Nq_WP; idx_q_bin++){

                    /* Upper and lower boundaries of current q-bin */
                    double q_bin_lower = q_WP_array[2*idx_q_bin    ];
                    double q_bin_upper = q_WP_array[2*idx_q_bin + 1];

                    cdouble R = {0, 0};
                    cdouble Q = {0, 0};
                    if (bound_state_exists){    // Calculate bound-continuum (BC) resolvent part R
                        R = resolvent_bound_continuum(E,
                                                      q_bin_upper,
                                                      q_bin_lower,
                                                      Eb);
                    }
                    else{                       // Calculate continuum-continuum (CC) resolvent part Q
                        Q = resolvent_continuum_continuum(E,
                                                          q_bin_upper,
                                                          q_bin_lower,
                                                          p_bin_upper,
                                                          p_bin_lower);
                    }

                    /* Use identical indexing as used in permutation matrix */
                    int G_idx = (idx_alpha-alpha_idx_lower)*Nq_WP*Np_WP + idx_q_bin*Np_WP + idx_p_bin;
                    G_subarray[G_idx] = R + Q;
                }
            }
        }
    }
}
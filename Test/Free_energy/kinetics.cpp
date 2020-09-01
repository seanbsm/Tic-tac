
#include "kinetics.h"

double calculate_3N_kinetic_energy(double* state_3N_bra_array,
                                   double* state_3N_ket_array,
                                   int &Np, double* p_array, double* wp_array,
                                   int &Nq, double* q_array, double* wq_array,
                                   int& Nalpha){

    double p, q, wp, wq, psi_bra, psi_ket, T;

    /* Kinetic energy test */
    double kinetic_energy_unormalised = 0;
    double inner_product  = 0;
    for (int idx_p=0; idx_p<Np; idx_p++){
        p  = p_array[idx_p];
        wp = wp_array[idx_p];
        
        for (int idx_q=0; idx_q<Nq; idx_q++){
            q  = q_array[idx_q];
            wq = wq_array[idx_q];
            
            T = (1./MN)*(p*p + 0.75*q*q);

            for (int idx_alpha=0; idx_alpha<Nalpha; idx_alpha++){
                //psi_symm = state_3N_symm_array[idx_p*Nq*Nalpha + idx_q * Nalpha + idx_alpha];
                psi_bra = state_3N_bra_array[idx_p*Nq*Nalpha + idx_q * Nalpha + idx_alpha];
                psi_ket = state_3N_ket_array[idx_p*Nq*Nalpha + idx_q * Nalpha + idx_alpha];
                
                kinetic_energy_unormalised +=  p*p*wp * q*q*wq  * T * psi_bra*psi_ket;
                inner_product += p*p*wp * q*q*wq * psi_bra*psi_ket;
            }
        }
    }
	
    return hbarc*hbarc*kinetic_energy_unormalised/inner_product;
}

double calculate_3N_potential_energy(double* state_3N_bra_array,
                                     double* state_3N_ket_array,
                                     int &Np, double* p_array, double* wp_array,
                                     int &Nq, double* q_array, double* wq_array,
                                     int& Nalpha){
	
    return 0;
}
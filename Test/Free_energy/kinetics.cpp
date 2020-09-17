
#include "kinetics.h"

double calculate_3N_kinetic_energy(double* state_3N_bra_array,
                                   double* state_3N_ket_array,
                                   int &Np, double* p_array, double* wp_array,
                                   int &Nq, double* q_array, double* wq_array,
                                   int& Nalpha, int* L_2N, int* S_2N, int* J_2N, int* T_2N, int* l_3N, int* two_j_3N){

    double p, q, wp, wq, psi_bra, psi_ket, kin_term;

    /* Kinetic energy test */
    double kinetic_energy_unormalised = 0;
    double inner_product  = 0;
    for (int idx_p=0; idx_p<Np; idx_p++){
        p  = p_array[idx_p];
        wp = wp_array[idx_p];
        
        for (int idx_q=0; idx_q<Nq; idx_q++){
            q  = q_array[idx_q];
            wq = wq_array[idx_q];
            
            kin_term = (1./MN)*(p*p + 0.75*q*q);

            for (int idx_alpha=0; idx_alpha<Nalpha; idx_alpha++){
                
                psi_bra = state_3N_bra_array[idx_p*Nq*Nalpha + idx_q * Nalpha + idx_alpha];
                psi_ket = state_3N_ket_array[idx_p*Nq*Nalpha + idx_q * Nalpha + idx_alpha];
                
                kinetic_energy_unormalised +=  p*p*wp * q*q*wq  * kin_term * psi_bra*psi_ket;
                inner_product += p*p*wp * q*q*wq * psi_bra*psi_ket;
            }
        }
    }
	
    return hbarc*hbarc*kinetic_energy_unormalised/inner_product;
}

double extract_potential_element_from_array(int& L, int& Lp, int& J, int& S, double* V_array){

    if (J<abs(S-L) or J>S+L){exit(1);}
    if (J<abs(S-Lp) or J>S+Lp){exit(1);}

    double potential_element = NAN;

    bool coupled = ( ((L!=Lp) or (L==Lp and L!=J)) );

    if (coupled){
        if (L==Lp and L<J){         // --
            potential_element = V_array[2];
        }
        else if (L==Lp and L>J){    // ++
            potential_element = V_array[5];
        }
        else if (L<Lp){             // +-
            potential_element = (-1)*V_array[3];
        }
        else{                       // -+
            potential_element = (-1)*V_array[4];
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

double calculate_3N_potential_energy(double* state_3N_bra_array,
                                     double* state_3N_ket_array,
                                     int &Np, double* p_array, double* wp_array,
                                     int &Nq, double* q_array, double* wq_array,
                                     int& Nalpha, int* L_2N, int* S_2N, int* J_2N, int* T_2N, int* l_3N, int* two_j_3N,
                                     potential_model* pot_ptr_np,potential_model* pot_ptr_nn){
    
    double p, pp, q, wp, wpp, wq, psi_bra, psi_ket, pot_term, pi, po;
    double V_nn_array [6];
    double V_np_array [6];

    int L, S, J, T, l3, tj3;
    int Lp, Sp, Jp, Tp, l3p, tj3p;
    bool coupled;

    /* Kinetic energy test */
    double potential_energy_unormalised = 0;
    double inner_product  = 0;

    for (int idx_alpha=0; idx_alpha<Nalpha; idx_alpha++){
        L   = L_2N[idx_alpha];
        S   = S_2N[idx_alpha];
        J   = J_2N[idx_alpha];
        T   = T_2N[idx_alpha];
        l3  = l_3N[idx_alpha];
        tj3 = two_j_3N[idx_alpha];

        std::cout << idx_alpha << " / " << Nalpha << std::endl;

        for (int idx_alpha_p=0; idx_alpha_p<Nalpha; idx_alpha_p++){
            Lp   = L_2N[idx_alpha_p];
            Sp   = S_2N[idx_alpha_p];
            Jp   = J_2N[idx_alpha_p];
            Tp   = T_2N[idx_alpha_p];
            l3p  = l_3N[idx_alpha_p];
            tj3p = two_j_3N[idx_alpha_p];

            /* Make sure states can couple via interaction */
            if (S==Sp and J==Jp and T==Tp and (abs(L-Lp)==2 or abs(L-Lp)==0) and l3==l3p and tj3==tj3p){
                
                coupled = ( (L!=Lp) or (L==Lp and L!=J) );

                for (int idx_p=0; idx_p<Np; idx_p++){
                    p  = p_array[idx_p];
                    wp = wp_array[idx_p];

                    pi = p*hbarc;
                    
                    for (int idx_pp=0; idx_pp<Np; idx_pp++){
                        pp  = p_array[idx_pp];
                        wpp = wp_array[idx_pp];

                        po = pp*hbarc;
                        
                        if (T==1){
                            pot_ptr_nn->V(pi, po, coupled, S, J, T, V_nn_array);
                            pot_ptr_np->V(pi, po, coupled, S, J, T, V_np_array);

                            pot_term = (1./3)*extract_potential_element_from_array(L, Lp, J, S, V_np_array)
                                     + (2./3)*extract_potential_element_from_array(L, Lp, J, S, V_nn_array);
                        }
                        else{
                            pot_ptr_np->V(pi, po, coupled, S, J, T, V_np_array);

                            pot_term = extract_potential_element_from_array(L, Lp, J, S, V_np_array);
                        }

                        for (int idx_q=0; idx_q<Nq; idx_q++){
                            q  = q_array[idx_q];
                            wq = wq_array[idx_q];

                            psi_ket = state_3N_ket_array[idx_p*Nq*Nalpha + idx_q * Nalpha + idx_alpha];
                            psi_bra = state_3N_bra_array[idx_pp*Nq*Nalpha + idx_q * Nalpha + idx_alpha_p];

                            potential_energy_unormalised += p*p*wp * pp*pp*wpp * q*q*wq * pot_term * psi_bra*psi_ket;
                        }
                    }
                }
            }
        }
    }

    
    for (int idx_p=0; idx_p<Np; idx_p++){
        p  = p_array[idx_p];
        wp = wp_array[idx_p];
        
        for (int idx_q=0; idx_q<Nq; idx_q++){
            q  = q_array[idx_q];
            wq = wq_array[idx_q];

            for (int idx_alpha=0; idx_alpha<Nalpha; idx_alpha++){
                psi_bra = state_3N_bra_array[idx_p*Nq*Nalpha + idx_q * Nalpha + idx_alpha];
                psi_ket = state_3N_ket_array[idx_p*Nq*Nalpha + idx_q * Nalpha + idx_alpha];
                
                inner_product += p*p*wp * q*q*wq * psi_bra*psi_ket;
            }
        }
    }

    std::cout << inner_product << std::endl;

    const double asymm_factor = 3;

    return hbarc*hbarc*hbarc*asymm_factor*potential_energy_unormalised / (inner_product);
}

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
                                     int& Nalpha, int* L_2N, int* S_2N, int* J_2N, int* T_2N, int* l_3N, int* two_j_3N,
                                     potential_model* pot_ptr_np,potential_model* pot_ptr_nn){
    
    double p, pp, q, wp, wpp, wq, psi_bra, psi_ket, pot_term, CG_coeff;
    double V_array [6];

    int L, S, J, T;
    int Lp, Sp, Jp, Tp;
    bool coupled;

    /* Triton ground state isospin */
    double T3 = 1./2, T3z = -1./2;

    /* Kinetic energy test */
    double potential_energy_unormalised = 0;
    double inner_product  = 0;

    for (int idx_alpha=0; idx_alpha<Nalpha; idx_alpha++){
        L = L_2N[idx_alpha];
        S = S_2N[idx_alpha];
        J = J_2N[idx_alpha];
        T = T_2N[idx_alpha];

        std::cout << idx_alpha << " / " << Nalpha << std::endl;

        for (int idx_alpha_p=0; idx_alpha_p<Nalpha; idx_alpha_p++){
            Lp = L_2N[idx_alpha_p];
            Sp = S_2N[idx_alpha_p];
            Jp = J_2N[idx_alpha_p];
            Tp = T_2N[idx_alpha_p];

            /* Make sure states can couple via interaction */
            if (S==Sp and J==Jp and T==Tp and abs(L-Lp)<3){
                coupled = ( ((L!=Lp) or (L==Lp and L!=J)) and J!=0);

                for (int idx_p=0; idx_p<Np; idx_p++){
                    p  = p_array[idx_p];
                    wp = wp_array[idx_p];
                    
                    for (int idx_pp=0; idx_pp<Np; idx_pp++){
                        pp  = p_array[idx_pp];
                        wpp = wp_array[idx_pp];

                        pot_term = 0;
                        /* Expand into isospin m-scheme from isospin j-scheme */
                        for (int T2z=-T; T2z<=T; T2z++){

                            /* Call correct potential class depending on nn-, np-, pp-interaction */
                            /* CGcoeff(double J, double m, double J1, double m1, double J2, double m2); */
                            if (T2z==-1){       // nn-pair, lone p-particle
                                pot_ptr_nn->V(p, pp, coupled, S, J, T, V_array);
                                CG_coeff = CGcoeff(T3, T3z, T, T2z, 1./2, +1./2);
                            }
                            else if (T2z==0){   // np-pair, lone n-particle
                                pot_ptr_np->V(p, pp, coupled, S, J, T, V_array);

                               // if (S==0 and J==0 and L==0 and idx_pp==0){
                               //     std::cout << p << " " << V_array[0] << std::endl;
                               // }

                                CG_coeff = CGcoeff(T3, T3z, T, T2z, 1./2, -1./2);
                            }
                            else{ /* We ignore pp-interactions for now (i.e. only simulate systems with a single proton) */
                                continue;
                            }
                            // ss 00 mm pm mp pp
                            /* Figure out which coupling it is */
                            
                            if (coupled){
                                if (L==Lp and L<J){
                                    pot_term += CG_coeff*CG_coeff*V_array[2];
                                    //std::cout << "2: " << V_array[2] << std::endl;
                                }
                                else if (L==Lp and L>J){
                                    pot_term += CG_coeff*CG_coeff*V_array[5];
                                    //std::cout << "5: " << V_array[5] << std::endl;
                                }
                                else if (L<Lp){
                                    pot_term += CG_coeff*CG_coeff*V_array[4];
                                    //std::cout << "4: " << V_array[4] << std::endl;
                                }
                                else{
                                    pot_term += CG_coeff*CG_coeff*V_array[3];
                                    //std::cout << "3: " << V_array[3] << std::endl;
                                }
                            }
                            else{
                                if (J==0 and L!=J){
                                    pot_term += CG_coeff*CG_coeff*V_array[5];
                                    //std::cout << "5: " << V_array[5] << std::endl;
                                }
                                else if (S==0){
                                    pot_term += CG_coeff*CG_coeff*V_array[0];
                                    //std::cout << "0: " << V_array[0] << std::endl;
                                }
                                else{
                                    pot_term += CG_coeff*CG_coeff*V_array[1];
                                    //std::cout << "1: " << V_array[1] << std::endl;
                                }
                            }
                        }

                        for (int idx_q=0; idx_q<Nq; idx_q++){
                            q  = q_array[idx_q];
                            wq = wq_array[idx_q];

                            psi_ket = state_3N_ket_array[idx_p*Nq*Nalpha + idx_q * Nalpha + idx_alpha];
                            psi_bra = state_3N_bra_array[idx_pp*Nq*Nalpha + idx_q * Nalpha + idx_alpha_p];

                            potential_energy_unormalised +=  p*p*wp * pp*pp*wpp * q*q*wq * pot_term * psi_bra*psi_ket;
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

    return hbarc*hbarc*M_PI*MN*0.5*potential_energy_unormalised/ (inner_product);
}
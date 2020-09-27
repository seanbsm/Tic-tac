
#include "state_antisymmetrization.h"

void antisymmetrize_state(double* state_3N_symm_array,
                          double* state_3N_asym_array,
                          double* A123,
                          int Np, int Nq, int Nalpha){
    
    //double state_symm, state_asym;

    int D123_dim = Np * Nq * Nalpha;

    //int A123_idx;
    
    #pragma omp parallel for 
    for (int alpha_idx=0; alpha_idx<Nalpha; alpha_idx++){
        double state_symm, state_asym;
        int A123_idx;
        for (int q_idx=0; q_idx<Nq; q_idx++){
            for (int p_idx=0; p_idx<Np; p_idx++){

                state_asym = 0;
                for (int alphap_idx=0; alphap_idx<Nalpha; alphap_idx++){
                    for (int qp_idx=0; qp_idx<Nq; qp_idx++){
                        for (int pp_idx=0; pp_idx<Np; pp_idx++){
                            
                            A123_idx = (int) (alpha_idx*Nq*Np + q_idx*Np + p_idx)*D123_dim +
                                               alphap_idx*Nq*Np + qp_idx*Np + pp_idx;

                            //A123_idx = (int) (alphap_idx*Nq*Np + qp_idx*Np + pp_idx)*D123_dim +
                            //                  alpha_idx*Nq*Np + q_idx*Np + p_idx;

                            state_symm = state_3N_symm_array[pp_idx*Nq*Nalpha + qp_idx * Nalpha + alphap_idx];

                            /* Essentially just a matrix-vector multiplication: psi = A123 x psi' */
                            state_asym += A123[A123_idx] * state_symm;
                        }
                    }
                }

                state_3N_asym_array[p_idx*Nq*Nalpha + q_idx*Nalpha + alpha_idx] = state_asym;
            }
        }
    }
}
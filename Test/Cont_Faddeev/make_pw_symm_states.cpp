
#include "make_pw_symm_states.h"

void construct_symmetric_pw_states(int two_J_3N,
                                   int two_T_3N,
                                   int parity_3N,
                                   int two_J_1N_min,
                                   int two_J_1N_max,
                                   int J_2N_min,
                                   int J_2N_max,
                                   int& Nalpha,
                                   int** L_2N_array_ptr,
                                   int** S_2N_array_ptr,
                                   int** J_2N_array_ptr,
                                   int** T_2N_array_ptr,
                                   int** l_3N_array_ptr,
                                   int** two_j_3N_array_ptr){

    bool print_content = false;

    if (print_content){
        /* Print partial-wave expansion (PWE) truncations */
        std::cout << "Truncating partial-wave (pw) state space with: \n"
                  << "two_J_3N:     " << two_J_3N      << "\n"
                  << "two_T_3N:     " << two_T_3N      << "\n"
                  << "parity_3N:    " << parity_3N     << "\n"
                  << "two_J_1N_min: " << two_J_1N_min  << "\n"
                  << "two_J_1N_max: " << two_J_1N_max  << "\n"
                  << "J_2N_min:     " << J_2N_min      << "\n"
                  << "J_2N_max:     " << J_2N_max      << "\n"
                  << std::endl;
    }

    /* Make sure the pw state space truncations are positive */
    if (two_J_1N_min<0 or two_J_1N_max<0){
        raise_error("Cannot have negative two_J_1N!");
    }
    if (J_2N_min<0 or J_2N_max<0){
        raise_error("Cannot have negative J_2N!");
    }
    
    /* Define and set parity validation parameter from input */
    int parity_exponential_remainder = 0;
    if (parity_3N>0){
        parity_exponential_remainder = 0;
    }
    else{
        parity_exponential_remainder = 1;
    }

    std::vector<int> L_2N_temp;
    std::vector<int> S_2N_temp;
    std::vector<int> J_2N_temp;
    std::vector<int> T_2N_temp;
    std::vector<int> l_3N_temp;
    std::vector<int> two_j_3N_temp;

    int Nalpha_temp = 0;
    int L_1N_min = 0;
    int L_1N_max = 0;
    int L_2N_min = 0;
    int L_2N_max = 0;
    /* J_2N loop */
    for (int J_2N=J_2N_min; J_2N<J_2N_max+1; J_2N++){
        /* S_2N loop */
        for (int S_2N=0; S_2N<2; S_2N++){
            L_2N_min = J_2N - S_2N;
            L_2N_max = J_2N + S_2N;
            /* L_2N loop */
            for (int L_2N=L_2N_min; L_2N<L_2N_max+1; L_2N++){
                /* L_2N can't be negative */
                if (L_2N >= 0){
                    /* T_2N loop */
                    for (int T_2N=0; T_2N<2; T_2N++){
                        /* Generalised Pauli exclusion principle: L_2N + S_2N + T_2N must be odd */
                        if ( (L_2N + S_2N + T_2N)%2 ){
                            for (int two_J_1N=two_J_1N_min; two_J_1N<two_J_1N_max+1; two_J_1N++){
                                /* Cannot have even two_J_1N (l is integer (i.e. 2n/2), s is 1/2 -> 2j = 2n+1) */
                                if ( two_J_1N%2 ){
                                    /* Check if triangle inequality is satisfied for two_J_3N */
                                    if ( std::abs(2*J_2N-two_J_1N)<=two_J_3N and two_J_3N<=(2*J_2N+two_J_1N) ){
                                        /* Check if triangle inequality is satisfied for two_T_3N */
                                        if ( std::abs(2*T_2N-1)<=two_T_3N and two_T_3N<=(2*T_2N+1) ){
                                            L_1N_min = (int) (two_J_1N - 1)/2;
                                            L_1N_max = (int) (two_J_1N + 1)/2;
                                            /* L_1N loop */
                                            for (int L_1N=L_1N_min; L_1N<L_1N_max+1; L_1N++){
                                                /* L_2N can't be negative */
                                                if (L_1N >= 0){
                                                    /* Check 3N-system total parity given by parity_3N */
                                                    if ( ((L_2N+L_1N)%2)==parity_exponential_remainder ){

                                                        /* We've found a physical state.
                                                         * Append to temporary vectors */
                                                        L_2N_temp.push_back(L_2N);
                                                        S_2N_temp.push_back(S_2N);
                                                        J_2N_temp.push_back(J_2N);
                                                        T_2N_temp.push_back(T_2N);
                                                        l_3N_temp.push_back(L_1N);
                                                        two_j_3N_temp.push_back(two_J_1N);
                                                        
                                                        if (print_content){
                                                            std::cout << Nalpha_temp << " "
                                                                      << L_2N << " "
                                                                      << S_2N << " "
                                                                      << J_2N << " "
                                                                      << T_2N << " "
                                                                      << L_1N << " "
                                                                      << two_J_1N << std::endl;
                                                        }

                                                        /* Increment state counter */
                                                        Nalpha_temp += 1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    /* Write number of states found to input integer */
    Nalpha = Nalpha_temp;

    /* Allocate arrays to input array pointers */
    *L_2N_array_ptr     = new int [Nalpha];
    *S_2N_array_ptr     = new int [Nalpha];
    *J_2N_array_ptr     = new int [Nalpha];
    *T_2N_array_ptr     = new int [Nalpha];
    *l_3N_array_ptr     = new int [Nalpha];
    *two_j_3N_array_ptr = new int [Nalpha];

    /* Write temporary vector contents to newly allocated arrays */
    std::copy( L_2N_temp.begin(), L_2N_temp.end(), *L_2N_array_ptr );
    std::copy( S_2N_temp.begin(), S_2N_temp.end(), *S_2N_array_ptr );
    std::copy( J_2N_temp.begin(), J_2N_temp.end(), *J_2N_array_ptr );
    std::copy( T_2N_temp.begin(), T_2N_temp.end(), *T_2N_array_ptr );
    std::copy( l_3N_temp.begin(), l_3N_temp.end(), *l_3N_array_ptr );
    std::copy( two_j_3N_temp.begin(), two_j_3N_temp.end(), *two_j_3N_array_ptr );
}
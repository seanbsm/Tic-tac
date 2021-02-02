
#include "make_pw_symm_states.h"

void construct_symmetric_pw_states(int parity_3N,
                                   int J_2N_max,
                                   int two_J_3N_max,
                                   int& Nalpha,
                                   int** L_2N_array_ptr,
                                   int** S_2N_array_ptr,
                                   int** J_2N_array_ptr,
                                   int** T_2N_array_ptr,
                                   int** L_1N_array_ptr,
                                   int** two_J_1N_array_ptr,
                                   int** two_J_3N_array_ptr,
                                   int** two_T_3N_array_ptr){

    bool print_content = true;
    char print_table_format_words[] = "%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s\n";
    char print_table_format_ints[]  = "%-10d%-10d%-10d%-10d%-10d%d/%-8d%-10d%d/%-8d%d/%-8d\n";

    if (print_content){
        /* Print partial-wave expansion (PWE) truncations */
        std::cout << "Truncating partial-wave (pw) state space with: \n"
                  << "parity_3N:    " << parity_3N     << "\n"
                  << "J_2N_max:     " << J_2N_max      << "\n"
                  << "two_J_3N_max: " << two_J_3N_max  << "\n"
                  << std::endl;
        
        printf(print_table_format_words, "alpha", "L_2N", "S_2N", "J_2N", "L_1N", "J_1N", "T_2N", "T_3N", "J_3N");
    }

    /* Make sure the pw state space truncations are positive */
    if (J_2N_max<0){
        raise_error("Cannot have negative J_2N!");
    }
    
    /* Define and set parity validation parameter from input */
    int parity_exponential_remainder = parity_3N%2;

    std::vector<int> L_2N_temp;
    std::vector<int> S_2N_temp;
    std::vector<int> J_2N_temp;
    std::vector<int> T_2N_temp;
    std::vector<int> L_1N_temp;
    std::vector<int> two_J_1N_temp;
    std::vector<int> two_J_3N_temp;
    std::vector<int> two_T_3N_temp;

    int Nalpha_temp = 0;
    int L_1N_min = 0;
    int L_1N_max = 0;
    int L_2N_min = 0;
    int L_2N_max = 0;
    int T_2N_min = 0;
    int T_2N_max = 0;
    int two_J_1N_min = 0;
    int two_J_1N_max = 0;
    /* two_J_3N loop */
    for (int two_J_3N=1; two_J_3N<two_J_3N_max+1; two_J_3N+=2){
        /* two_T_3N loop */
        for (int two_T_3N=1; two_T_3N<4; two_T_3N+=2){
            /* J_2N loop */
            for (int J_2N=0; J_2N<J_2N_max+1; J_2N++){
                /* S_2N loop */
                for (int S_2N=0; S_2N<2; S_2N++){
                    /* Triangle inequality for L_2N */
                    L_2N_min = (int) abs(J_2N - S_2N);
                    L_2N_max = J_2N + S_2N;
                    /* L_2N loop */
                    for (int L_2N=L_2N_min; L_2N<L_2N_max+1; L_2N++){
                        /* Triangle inequality for T_2N */
                        T_2N_min = (int) abs(two_T_3N - 1)/2;
                        T_2N_max = (two_T_3N + 1)/2;
                        /* T can't be greater than 1 */
                        if (T_2N_max>1){
                            T_2N_max=1;
                        }
                        /* T_2N loop */
                        for (int T_2N=T_2N_min; T_2N<T_2N_max+1; T_2N++){
                            /* Generalised Pauli exclusion principle: L_2N + S_2N + T_2N must be odd */
                            if ( (L_2N + S_2N + T_2N)%2 ){
                                /* Triangle inequality for J_1N */
                                two_J_1N_min = (int) abs(two_J_3N - 2*J_2N);
                                two_J_1N_max = two_J_3N + 2*J_2N;
                                /* two_J_1N loop */
                                for (int two_J_1N=two_J_1N_min; two_J_1N<two_J_1N_max+1; two_J_1N+=2){
                                    /* Triangle inequality for L_1N */
                                    L_1N_min = (int) (abs(two_J_1N - 1))/2;
                                    L_1N_max = (int) (two_J_1N + 1)/2;
                                    /* L_1N loop */
                                    for (int L_1N=L_1N_min; L_1N<L_1N_max+1; L_1N++){
                                        /* Check 3N-system total parity given by parity_3N */
                                        //if ( ((L_2N+L_1N)%2)==parity_exponential_remainder ){
                                            
                                            //if (two_T_3N==3){
                                            //    if ( (S_2N==0 && L_2N==0 && J_2N==0)==false ){
                                            //        continue;
                                            //    }
                                            //}

                                            /* We've found a physical state.
                                             * Append to temporary vectors */
                                            L_2N_temp.push_back(L_2N);
                                            S_2N_temp.push_back(S_2N);
                                            J_2N_temp.push_back(J_2N);
                                            T_2N_temp.push_back(T_2N);
                                            L_1N_temp.push_back(L_1N);
                                            two_J_1N_temp.push_back(two_J_1N);
                                            two_J_3N_temp.push_back(two_J_3N);
                                            two_T_3N_temp.push_back(two_T_3N);

                                            /* Prints in the same order as table 1 of Glockle et al., Phys. Rep. 274 (1996) 107-285 */
                                            if (print_content){
                                                printf (print_table_format_ints, Nalpha_temp, L_2N, S_2N, J_2N, L_1N, two_J_1N,2, T_2N, two_T_3N,2, two_J_3N,2);
                                            }

                                            /* Increment state counter */
                                            Nalpha_temp += 1;
                                        //}
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
    *L_1N_array_ptr     = new int [Nalpha];
    *two_J_1N_array_ptr = new int [Nalpha];
    *two_J_3N_array_ptr = new int [Nalpha];
    *two_T_3N_array_ptr = new int [Nalpha];

    /* Write temporary vector contents to newly allocated arrays */
    std::copy( L_2N_temp.begin(), L_2N_temp.end(), *L_2N_array_ptr );
    std::copy( S_2N_temp.begin(), S_2N_temp.end(), *S_2N_array_ptr );
    std::copy( J_2N_temp.begin(), J_2N_temp.end(), *J_2N_array_ptr );
    std::copy( T_2N_temp.begin(), T_2N_temp.end(), *T_2N_array_ptr );
    std::copy( L_1N_temp.begin(), L_1N_temp.end(), *L_1N_array_ptr );
    std::copy( two_J_1N_temp.begin(), two_J_1N_temp.end(), *two_J_1N_array_ptr );
    std::copy( two_J_3N_temp.begin(), two_J_3N_temp.end(), *two_J_3N_array_ptr );
    std::copy( two_T_3N_temp.begin(), two_T_3N_temp.end(), *two_T_3N_array_ptr );
}
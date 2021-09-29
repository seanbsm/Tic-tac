
#include "make_pw_symm_states.h"

int unique_2N_idx(int L_2N, int S_2N, int J_2N, int T_2N, bool coupled, run_params run_parameters){
	
	if (coupled==true && run_parameters.tensor_force==false){
		raise_error("Cannot have coupled states without a tensor force.");
	}

	/* If isospin-breaking in 1S0 is enables (isospin_breaking_1S0=true)
	 * we change the unique index to move 1S0 from an uncoupled to a coupled state.
	 * This is because with isospin-breaking 1S0 essentially becomes coupled via T_3N: 1/2 <-> 3/2 */
	
	bool state_1S0 = (S_2N==0 && J_2N==0 && L_2N==0 && T_2N==1);

	int unique_idx;
	if (run_parameters.tensor_force==true){
		if (coupled){
			/* Unique index for all coupled states if tensor force is on */
			unique_idx = J_2N-1;

			/* We give room to 1S0 as a coupled state if isospin-breaking is enabled.
			 * All other coupled states are moved up and 1S0 gets index 0. */
			if (state_1S0==false && run_parameters.isospin_breaking_1S0==true){
				unique_idx += 1;
			}
			if (state_1S0==true  && run_parameters.isospin_breaking_1S0==true){
				unique_idx  = 0;
			}
		}
		else{
			/* Unique index for all uncoupled 2N states if tensor force is on */
			unique_idx = 2*J_2N + S_2N;

			/* We remove the uncoupled slot given to 1S0 if it is coupled */
			if (run_parameters.isospin_breaking_1S0==true){
				unique_idx -= 1;
			}
		}
	}
	else{
		if (coupled){
			/* 1S0 will be the only coupled state if the tensor force is off */
			if (run_parameters.isospin_breaking_1S0==true){
				unique_idx = 0;
			}
		}
		else{
			/* Unique index for all uncoupled 2N states if tensor force is off */
			unique_idx = J_2N*(4*(J_2N>1) + 2) + (S_2N==1) + (S_2N==1 && J_2N!=0)*(J_2N-L_2N + 1);

			/* We remove the uncoupled slot given to 1S0 if it is coupled */
			if (run_parameters.isospin_breaking_1S0==true){
				unique_idx -= 1;
			}
		}
	}

	return unique_idx;
}

void construct_symmetric_pw_states(int&   N_chn_3N,
								   int**  chn_3N_idx_array_ptr,
								   int&   Nalpha,
								   int**  L_2N_array_ptr,
								   int**  S_2N_array_ptr,
								   int**  J_2N_array_ptr,
								   int**  T_2N_array_ptr,
								   int**  L_1N_array_ptr,
								   int**  two_J_1N_array_ptr,
								   int**  two_J_3N_array_ptr,
								   int**  two_T_3N_array_ptr,
								   int**  P_3N_array_ptr,
								   pw_3N_statespace& pw_states,
								   run_params run_parameters){
	
	int J_2N_max	 = run_parameters.J_2N_max;
	int two_J_3N_max = run_parameters.two_J_3N_max;
	
	int two_T_3N_max = 1;
	if (run_parameters.isospin_breaking_1S0==true){
		two_J_3N_max = 3;
	}

	bool print_content = true;
	char print_table_format_words[] = "%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s\n";
	char print_table_format_ints[]  = "%-10d%-10d%-10d%-10d%-10d%d/%-8d%-10d%d/%-8d%d/%-8d%-10d\n";

	//char* print_table_format_words = new char [80000];
	//char* print_table_format_ints  = new char [80000];
	//sprintf(print_table_format_words, "%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s\n");
	//sprintf(print_table_format_ints,  "%-10d%-10d%-10d%-10d%-10d%d/%-8d%-10d%d/%-8d%d/%-8d%-10d\n");

	if (print_content){
		/* Print partial-wave expansion (PWE) truncations */
		std::cout << "Truncating partial-wave (pw) state space with: \n"
				  << "J_2N_max:     " << J_2N_max      << "\n"
				  << "two_J_3N_max: " << two_J_3N_max  << "\n"
				  << std::endl;
		
		printf(print_table_format_words, "alpha", "L_2N", "S_2N", "J_2N", "L_1N", "J_1N", "T_2N", "T_3N", "J_3N", "PAR");
	}

	/* Make sure the pw state space truncations are positive */
	if (J_2N_max<0){
		raise_error("Cannot have negative J_2N!");
	}

	std::vector<int> chn_idx_temp;
	
	std::vector<int>  L_2N_temp;
	std::vector<int>  S_2N_temp;
	std::vector<int>  J_2N_temp;
	std::vector<int>  T_2N_temp;
	std::vector<int>  L_1N_temp;
	std::vector<int>  two_J_1N_temp;
	std::vector<int>  two_J_3N_temp;
	std::vector<int>  two_T_3N_temp;
	std::vector<int>  P_3N_temp;

	int  Nalpha_temp      = 0;
	int  Nalpha_temp_prev = 0;
	int  N_chn_3N_temp    = 0;
	int  L_1N_min         = 0;
	int  L_1N_max         = 0;
	int  L_2N_min         = 0;
	int  L_2N_max         = 0;
	int  T_2N_min         = 0;
	int  T_2N_max         = 0;
	int  two_J_1N_min     = 0;
	int  two_J_1N_max     = 0;
	int  P_3N             = 0;
	/* two_J_3N loop */
	for (int two_J_3N=1; two_J_3N<two_J_3N_max+1; two_J_3N+=2){
		/* Parity loop: P_3N_remainder is the remainder of (l+L)/2 */
		for (int P_3N_remainder=0; P_3N_remainder<2; P_3N_remainder++){
			if (print_content){
				if (P_3N_remainder==0){
					printf ("Constructing channel %d with JP=%d/2-: \n", N_chn_3N_temp+1, two_J_3N);
				}
				else{
					printf ("Constructing channel %d with JP=%d/2+: \n", N_chn_3N_temp+1, two_J_3N);
				}
			}
			/* Used to cound how many PW states are in current channel (J_3N, P_3N) */
			int Nalpha_in_current_chn = 0;
			/* two_T_3N loop */
			for (int two_T_3N=1; two_T_3N<two_J_3N_max+1; two_T_3N+=2){
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
											/* Check 3N-system total parity given by P_3N */
											if ( ((L_2N+L_1N)%2)==P_3N_remainder){

												if (two_T_3N==3){
												    if ( (S_2N==0 && L_2N==0 && J_2N==0)==false ){
												        continue;
												    }
												}

												/* Quartet channel for Malfliet-Tjon debugging */
												//if (L_2N!=0 ||
												//	S_2N!=0 ||
												//	J_2N!=0 ||
												//	T_2N!=1){
												//	continue;
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

												/* +1 if P_3N_remainder=0, -1 if P_3N_remainder=1 */
												P_3N = 1 - 2*P_3N_remainder;
												P_3N_temp.push_back(P_3N);

												/* Prints in the same order as table 1 of Glockle et al., Phys. Rep. 274 (1996) 107-285 */
												if (print_content){
													printf (print_table_format_ints, Nalpha_temp, L_2N, S_2N, J_2N, L_1N, two_J_1N,2, T_2N, two_T_3N,2, two_J_3N,2, P_3N);
												}

												/* Increment state counters */
												Nalpha_temp           += 1;
												Nalpha_in_current_chn += 1;
											}
										}
									}
								}
							}
						}
					}
				}
			}
			/* Append 3N channel if there exists states in the given channel (J_3N, T_3N, P_3N)
			 * Using Nalpha_temp_prev rather than Nalpha_temp means we get the starting index
			 * for the current channel. It makes for simpler indexing throughout the code. */
			if (Nalpha_in_current_chn!=0){
				N_chn_3N_temp += 1;
				chn_idx_temp.push_back(Nalpha_temp_prev);
			}
			Nalpha_temp_prev = Nalpha_temp;
		}
	}

	/* Append the state space size as well - this allows for simpler indexing */
	chn_idx_temp.push_back(Nalpha_temp);
	
	/* Write number of states found to input integer */
	Nalpha = Nalpha_temp;

	/* Write number of 3N channels found to input integer */
	N_chn_3N = N_chn_3N_temp;
	
	/* Allocate arrays to input array pointers */
	*L_2N_array_ptr     = new int [Nalpha];
	*S_2N_array_ptr     = new int [Nalpha];
	*J_2N_array_ptr     = new int [Nalpha];
	*T_2N_array_ptr     = new int [Nalpha];
	*L_1N_array_ptr     = new int [Nalpha];
	*two_J_1N_array_ptr = new int [Nalpha];
	*two_J_3N_array_ptr = new int [Nalpha];
	*two_T_3N_array_ptr = new int [Nalpha];
	*P_3N_array_ptr     = new int [Nalpha];

	*chn_3N_idx_array_ptr = new int [N_chn_3N+1];

	/* Write temporary vector contents to newly allocated arrays */
	std::copy( L_2N_temp.begin(), L_2N_temp.end(), *L_2N_array_ptr );
	std::copy( S_2N_temp.begin(), S_2N_temp.end(), *S_2N_array_ptr );
	std::copy( J_2N_temp.begin(), J_2N_temp.end(), *J_2N_array_ptr );
	std::copy( T_2N_temp.begin(), T_2N_temp.end(), *T_2N_array_ptr );
	std::copy( L_1N_temp.begin(), L_1N_temp.end(), *L_1N_array_ptr );
	std::copy( two_J_1N_temp.begin(), two_J_1N_temp.end(), *two_J_1N_array_ptr );
	std::copy( two_J_3N_temp.begin(), two_J_3N_temp.end(), *two_J_3N_array_ptr );
	std::copy( two_T_3N_temp.begin(), two_T_3N_temp.end(), *two_T_3N_array_ptr );
	std::copy( P_3N_temp.begin(), P_3N_temp.end(), *P_3N_array_ptr );

	std::copy( chn_idx_temp.begin(), chn_idx_temp.end(), *chn_3N_idx_array_ptr );

	pw_states.dim				=  Nalpha;
	pw_states.J_2N_max			=  J_2N_max;
	pw_states.L_2N_array		= *L_2N_array_ptr;
	pw_states.S_2N_array		= *S_2N_array_ptr;
	pw_states.J_2N_array		= *J_2N_array_ptr;
	pw_states.T_2N_array		= *T_2N_array_ptr;
	pw_states.L_1N_array		= *L_1N_array_ptr;
	pw_states.two_J_1N_array	= *two_J_1N_array_ptr;
	pw_states.two_J_3N_array	= *two_J_3N_array_ptr;
	pw_states.two_T_3N_array	= *two_T_3N_array_ptr;
	pw_states.P_3N_array		= *P_3N_array_ptr;
}
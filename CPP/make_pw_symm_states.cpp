
#include "make_pw_symm_states.h"

int make_2N_chn_idx(int L,  int S,  int J,  int T, int Jmax){
	int Lmax = Jmax+1;
	int Smax = 1;
	int Tmax = 1;
	int num_2N_chns = Jmax*6 + 2;

	int idx  = J *(Lmax+1)*(Smax+1)*(Tmax+1) + L *(Smax+1)*(Tmax+1) + S *(Tmax+1) + T;
	
	return idx;
}

bool check_2N_coupling(int L,  int S,  int J,  int T,
					   int Lp, int Sp, int Jp, int Tp,
					   bool tensor_force_on){
	
	bool coupled = false;

	bool L_allowed = (L==Lp);
	if (tensor_force_on==true && L_allowed==false){
		L_allowed = abs(L-Lp)==2;
	}
	bool S_allowed = (S==Sp);
	bool J_allowed = (J==Jp);
	bool T_allowed = (T==Tp);

	if (L_allowed && S_allowed && J_allowed && T_allowed){
		coupled = true;
	}
	else{
		coupled = false;
	}

	return coupled;
}

void construct_symmetric_pw_states(int    J_2N_max,
								   int    two_J_3N_max,
								   int&   N_chn_3N,
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
								   bool   tensor_force_on,
								   int&   Nalpha_1N,
								   int&   Nalpha_2N,
								   pw_1N_state** pw_1N_states_array_ptr,
								   pw_2N_state** pw_2N_states_array_ptr,
								   pw_3N_state** pw_3N_states_array_ptr){

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
		/* two_T_3N loop */
		for (int two_T_3N=1; two_T_3N<2; two_T_3N+=2){
			/* Parity loop: P_3N_remainder is the remainder of (l+L)/2 */
			for (int P_3N_remainder=0; P_3N_remainder<2; P_3N_remainder++){
				/* Used to cound how many PW states are in current channel (J_3N, T_3N, P_3N) */
				int Nalpha_in_current_chn = 0;
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

												//if (two_T_3N==3){
												//    if ( (S_2N==0 && L_2N==0 && J_2N==0)==false ){
												//        continue;
												//    }
												//}

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


	/* THE CODE BELOW GENERATES STRUCTS FOR THE STATE SPACE CONSTRUCTED ABOVE */

	Nalpha_1N = 0;
	Nalpha_2N = 0;

	std::vector<pw_1N_state> pw_1N_states_temp;
	std::vector<pw_2N_state> pw_2N_states_temp;
	std::vector<pw_3N_state> pw_3N_states_temp;

	for (int idx_alpha=0; idx_alpha<Nalpha; idx_alpha++){
		/* Retrieve quantum numbers of 3N state */
		int L_2N     = L_2N_temp    [idx_alpha];
		int S_2N     = S_2N_temp    [idx_alpha];
		int J_2N     = J_2N_temp    [idx_alpha];
		int T_2N     = T_2N_temp    [idx_alpha];
		int L_1N     = L_1N_temp    [idx_alpha];
		int two_J_1N = two_J_1N_temp[idx_alpha];
		int two_J_3N = two_J_3N_temp[idx_alpha];
		int two_T_3N = two_T_3N_temp[idx_alpha];
		int P_3N	 = P_3N_temp	[idx_alpha];

		/* See if we've already appended 1N quantum numbers to 1N state space */
		bool state_1N_found = false;
		for (int idx_1N=0; idx_1N<Nalpha_1N; idx_1N++){
			int L_1N_p     = pw_1N_states_temp[idx_1N].L_1N;
			int two_J_1N_p = pw_1N_states_temp[idx_1N].two_J_1N;
			if (L_1N==L_1N_p &&
				two_J_1N==two_J_1N_p){
				state_1N_found = true;
				break;
			}
		}
		/* Append 1N state to 1N state space */
		if (state_1N_found==false){
			pw_1N_state state_1N_new;
			state_1N_new.L_1N     = L_1N;
			state_1N_new.two_J_1N = two_J_1N;
			state_1N_new.two_S_1N = 1;
			state_1N_new.two_T_1N = 1;
			state_1N_new.state_idx = Nalpha_1N;
			pw_1N_states_temp.push_back(state_1N_new);
			Nalpha_1N += 1;
		}

		/* See if we've already appended 2N quantum numbers to 2N state space */
		bool state_2N_found = false;
		for (int idx_2N=0; idx_2N<Nalpha_2N; idx_2N++){
			int L_2N_p = pw_2N_states_temp[idx_2N].L_2N;
			int J_2N_p = pw_2N_states_temp[idx_2N].J_2N;
			int S_2N_p = pw_2N_states_temp[idx_2N].S_2N;
			int T_2N_p = pw_2N_states_temp[idx_2N].T_2N;
			if (L_2N==L_2N_p &&
				J_2N==J_2N_p &&
				S_2N==S_2N_p &&
				T_2N==T_2N_p){
				state_2N_found = true;
				break;
			}
		}
		/* Append 2N state to 2N state space */
		if (state_2N_found==false){
			pw_2N_state state_2N_new;
			state_2N_new.L_2N = L_2N;
			state_2N_new.J_2N = J_2N;
			state_2N_new.S_2N = S_2N;
			state_2N_new.T_2N = T_2N;
			state_2N_new.state_idx = Nalpha_2N;
			pw_2N_states_temp.push_back(state_2N_new);
			Nalpha_2N += 1;
		}

		/* Append 3N state to 3N state space (one-to-one mapping -> no lookup needed) */
		pw_3N_state state_3N_new;
		state_3N_new.L_1N	   = L_1N;
		state_3N_new.two_S_1N  = 1;
		state_3N_new.two_J_1N  = two_J_1N;
		state_3N_new.two_T_1N  = 1;
		state_3N_new.L_2N	   = L_2N;
		state_3N_new.S_2N	   = S_2N;
		state_3N_new.J_2N	   = J_2N;
		state_3N_new.T_2N	   = T_2N;
		state_3N_new.two_J_3N  = two_J_3N;
		state_3N_new.two_T_3N  = two_T_3N;
		state_3N_new.P_3N	   = P_3N;
		state_3N_new.state_idx = idx_alpha;
		pw_3N_states_temp.push_back(state_3N_new);
	}

	*pw_1N_states_array_ptr = new pw_1N_state [Nalpha_1N];
	*pw_2N_states_array_ptr = new pw_2N_state [Nalpha_2N];
	*pw_3N_states_array_ptr = new pw_3N_state [Nalpha];

	std::copy( pw_1N_states_temp.begin(), pw_1N_states_temp.end(), *pw_1N_states_array_ptr );
	std::copy( pw_2N_states_temp.begin(), pw_2N_states_temp.end(), *pw_2N_states_array_ptr );
	std::copy( pw_3N_states_temp.begin(), pw_3N_states_temp.end(), *pw_3N_states_array_ptr );

	/* Assign pointers to 1N and 2N states for a given 3N state */
	for (int idx_alpha_3N=0; idx_alpha_3N<Nalpha; idx_alpha_3N++){
		pw_3N_state state_3N = (*pw_3N_states_array_ptr)[idx_alpha_3N];

		for (int idx_alpha_1N=0; idx_alpha_1N<Nalpha_1N; idx_alpha_1N++){
			pw_1N_state state_1N = *pw_1N_states_array_ptr[idx_alpha_1N];
			std::cout << "L" << std::endl;
			std::cout << state_1N.L_1N << std::endl;
			if (state_1N.L_1N == state_3N.L_1N &&
				state_1N.two_J_1N == state_3N.two_J_1N){
				(*pw_3N_states_array_ptr)[idx_alpha_3N].pw_1N_state_idx = idx_alpha_1N;
			}
		}

		for (int idx_alpha_2N=0; idx_alpha_2N<Nalpha_2N; idx_alpha_2N++){
			pw_2N_state state_2N = (*pw_2N_states_array_ptr)[idx_alpha_2N];
			if (state_2N.L_2N == state_3N.L_2N &&
				state_2N.S_2N == state_3N.S_2N &&
				state_2N.J_2N == state_3N.J_2N &&
				state_2N.T_2N == state_3N.T_2N){
				(*pw_3N_states_array_ptr)[idx_alpha_3N].pw_2N_state_idx = idx_alpha_2N;
			}
		}
	}

	/* Search for coupling of 2N states */
	for (int idx_alpha_2N_r=0; idx_alpha_2N_r<Nalpha_2N; idx_alpha_2N_r++){
		pw_2N_state state_2N_r = (*pw_2N_states_array_ptr)[idx_alpha_2N_r];

		/* Initially set the coupling condition to false */
		state_2N_r.coupled 			 = false;
		state_2N_r.coupled_state_idx = idx_alpha_2N_r;

		/* Search for coupled states, given that tensor-force is allowed */
		if (tensor_force_on){

			/* Extract row-state quantum numbers */
			int L_2N_r = state_2N_r.L_2N;
			int S_2N_r = state_2N_r.S_2N;
			int J_2N_r = state_2N_r.J_2N;
			int T_2N_r = state_2N_r.T_2N;

			for (int idx_alpha_2N_c=0; idx_alpha_2N_c<Nalpha_2N; idx_alpha_2N_c++){
				pw_2N_state state_2N_c = (*pw_2N_states_array_ptr)[idx_alpha_2N_c];

				/* Extract row-state quantum numbers */
				int L_2N_c = state_2N_c.L_2N;
				int S_2N_c = state_2N_c.S_2N;
				int J_2N_c = state_2N_c.J_2N;
				int T_2N_c = state_2N_c.T_2N;

				if ( (abs(L_2N_r-L_2N_c)==0 || abs(L_2N_r-L_2N_c)==2) &&
					  S_2N_r==S_2N_c &&
					  J_2N_r==J_2N_c &&
					  T_2N_r==T_2N_c ){
					(*pw_2N_states_array_ptr)[idx_alpha_2N_r].coupled_state_idx = idx_alpha_2N_c;
					(*pw_2N_states_array_ptr)[idx_alpha_2N_c].coupled_state_idx = idx_alpha_2N_r;
				}
			}
		}
	}
}
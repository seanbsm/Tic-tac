#ifndef MAKE_PW_SYMM_STATES_H
#define MAKE_PW_SYMM_STATES_H

#include <iostream>
#include <vector>
#include <fstream>

#include "error_management.h"

struct pw_1N_state{
	int L_1N;
	int two_S_1N;
	int two_J_1N;
	int two_T_1N;
	int state_idx;
};

struct pw_2N_state{
	int L_2N;
	int S_2N;
	int J_2N;
	int T_2N;
	bool coupled;
	int coupled_state_idx;
	int state_idx;
};

struct pw_3N_state{
	int L_1N;
	int two_S_1N;
	int two_J_1N;
	int two_T_1N;
	int L_2N;
	int S_2N;
	int J_2N;
	int T_2N;
	int two_J_3N;
	int two_T_3N;
	int P_3N;
	int pw_1N_state_idx;
	int pw_2N_state_idx;
	int state_idx;
};

bool check_2N_coupling(int L,  int S,  int J,  int T,
					   int Lp, int Sp, int Jp, int Tp,
					   bool tensor_force_on);

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
								   pw_3N_state** pw_3N_states_array_ptr);

#endif // MAKE_PW_SYMM_STATES_H
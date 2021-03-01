#ifndef MAKE_PW_SYMM_STATES_H
#define MAKE_PW_SYMM_STATES_H

#include <iostream>
#include <vector>
#include <fstream>

#include "error_management.h"

void construct_symmetric_pw_states(int   J_2N_max,
								   int   two_J_3N_max,
								   int&  N_chn_3N,
								   int** chn_3N_idx_array_ptr,
								   int&  Nalpha,
								   int** L_2N_array_ptr,
								   int** S_2N_array_ptr,
								   int** J_2N_array_ptr,
								   int** T_2N_array_ptr,
								   int** L_1N_array_ptr,
								   int** two_J_1N_array_ptr,
								   int** two_J_3N_array_ptr,
								   int** two_T_3N_array_ptr,
								   int** P_3N_array_ptr);

#endif // MAKE_PW_SYMM_STATES_H
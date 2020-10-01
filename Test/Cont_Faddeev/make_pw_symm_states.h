#ifndef MAKE_PW_SYMM_STATES_H
#define MAKE_PW_SYMM_STATES_H

#include <iostream>
#include <vector>

#include "error_management.h"

void construct_symmetric_pw_states(int two_J_3N,
                                   int two_T_3N,
                                   int parity_3N,
                                   int& Nalpha,
                                   int** L_2N_array_ptr,
                                   int** S_2N_array_ptr,
                                   int** J_2N_array_ptr,
                                   int** T_2N_array_ptr,
                                   int** l_3N_array_ptr,
                                   int** two_j_3N_array_ptr);

#endif // MAKE_PW_SYMM_STATES_H
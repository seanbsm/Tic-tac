#ifndef MAKE_PW_SYMM_STATES_H
#define MAKE_PW_SYMM_STATES_H

#include <iostream>
#include <vector>
#include <fstream>

#include "error_management.h"
#include "type_defs.h"

int unique_2N_idx(int L_2N, int S_2N, int J_2N, int T_2N, bool coupled, run_params run_parameters);

void construct_symmetric_pw_states(pw_3N_statespace& pw_states,
								   run_params run_parameters);

#endif // MAKE_PW_SYMM_STATES_H
#ifndef STATE_ANTISYMMETRIZATION_H
#define STATE_ANTISYMMETRIZATION_H

#include <iostream>

void antisymmetrize_state(double* state_3N_symm_array,
                          double* state_3N_asym_array,
                          double* A123,
                          int Np, int Nq, int Nalpha);

#endif // STATE_ANTISYMMETRIZATION_H
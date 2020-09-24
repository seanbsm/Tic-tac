#ifndef PERMUTATION_OPERATORS_H
#define PERMUTATION_OPERATORS_H

#include <string>
#include "hdf5/serial/hdf5.h"
#include "hdf5/serial/hdf5_hl.h"

#include "Kai_code_and_data_package/read_psi_Sean.h"

void read_P123_data_file(double* P123);

void calculate_antisymmetrization_operator(double* state_3N_bra_array,
                                           double* state_3N_ket_array,
                                           int &Np, double* p_array, double* wp_array,
                                           int &Nq, double* q_array, double* wq_array,
                                           int& Nalpha, int* L_2N, int* S_2N, int* J_2N, int* T_2N, int* l_3N, int* two_j_3N,
                                           double* A123);

#endif // PERMUTATION_OPERATORS_H
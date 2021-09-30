#ifndef MAKE_POTENTIAL_MATRIX_H
#define MAKE_POTENTIAL_MATRIX_H

#include <iostream>
#include <math.h>

#include "error_management.h"

#include "type_defs.h"
#include "make_wp_states.h"
#include "make_pw_symm_states.h"
#include "Interactions/potential_model.h"

double extract_potential_element_from_array(int L, int Lp, int J, int S, bool coupled, double* V_array);

void calculate_potential_matrices_array_in_WP_basis(double* V_WP_unco_array, int num_2N_unco_states,
													double* V_WP_coup_array, int num_2N_coup_states,
													int Np_WP, double* p_WP_array,
													int Np, double* p_array, double* wp_array,
													int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array, int* two_T_3N_array,
													potential_model* pot_ptr_nn,
													potential_model* pot_ptr_np,
													run_params run_parameters);

#endif // MAKE_POTENTIAL_MATRIX_H
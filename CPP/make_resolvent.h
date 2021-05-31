#ifndef MAKE_RESOLVENT_H
#define MAKE_RESOLVENT_H

#include <iostream>

#include "error_management.h"
#include "constants.h"
#include "make_pw_symm_states.h"
#include "General_functions/kinetic_conversion.h"

/* COMMENTARY A:
 * Taken from Kukulin 2007 (https://doi.org/10.1007/s11232-007-0030-3), Eqs. (32), (34-35).
 * In their notation (left: program variables, right: LaTeX symbols from paper):
 * Eq_upper = \Epsilon_{j}
 * Eq_lower = \Epsilon_{j-1}
 * Ep_upper = E_{i}
 * Ep_lower = E_{i-1}
 * Dq       = \Delta_{j}
 * Dp       = \delta_{i}
 * DP       = \Delta_{+}
 * DM       = \Delta_{-}
 * D        = \Delta */

double heaviside_step_function(double val);

cdouble resolvent_bound_continuum(double E, double Eb,
								  double q_bin_upper,
								  double q_bin_lower);

cdouble resolvent_continuum_continuum(double E, double Eb,
									  double q_bin_upper,
									  double q_bin_lower,
									  double e_bin_upper,
									  double e_bin_lower);

void calculate_resolvent_array_in_SWP_basis(cdouble* G_array,
											double   E,
											bool 	 tensor_force_true,
											int      Np_WP, double* e_SWP_unco_array, double* e_SWP_coup_array,
											int      Nq_WP, double* q_WP_array,
											int      Nalpha,
											int*     L_2N_array,
											int*     S_2N_array,
											int*     J_2N_array,
											int*     T_2N_array);

#endif // MAKE_RESOLVENT_H
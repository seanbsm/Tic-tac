#ifndef MAKE_RESOLVENT_H
#define MAKE_RESOLVENT_H

#include <iostream>

#include "error_management.h"
#include "constants.h"

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

cdouble resolvent_bound_continuum(double E,
								  double q_bin_upper,
								  double q_bin_lower,
								  double Eb);

cdouble resolvent_continuum_continuum(double E,
									  double q_bin_upper,
									  double q_bin_lower,
									  double p_bin_upper,
									  double p_bin_lower);

void calculate_resolvent_array_in_SWP_basis(cdouble* G_array,
											double   E,
											int      Np_WP, double* p_SWP_unco_array, double* p_SWP_coup_array,
											int      Nq_WP, double* q_WP_array,
											int      Nalpha,
											int*     L_2N_array,
											int*     S_2N_array,
											int*     J_2N_array,
											int*     T_2N_array);

#endif // MAKE_RESOLVENT_H
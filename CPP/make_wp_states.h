#ifndef MAKE_WP_STATES_H
#define MAKE_WP_STATES_H

#include <iostream>
#include <math.h>

#include "constants.h"
#include "error_management.h"
#include "disk_io_routines.h"
#include "General_functions/kinetic_conversion.h"
#include "General_functions/gauss_legendre.h"

double p_normalization(double p0, double p1);
double q_normalization(double q0, double q1);

double p_weight_function(double p);
double q_weight_function(double q);

void make_chebyshev_distribution(int N_WP,
								 double* boundary_array,
								 double scale,
								 double	sparseness_degree);

void make_p_bin_grid(fwp_statespace& fwp_states, run_params run_parameters);
void make_q_bin_grid(fwp_statespace& fwp_states, run_params run_parameters);

void make_p_bin_quadrature_grids(fwp_statespace& fwp_states);
void make_q_bin_quadrature_grids(fwp_statespace& fwp_states);

void make_fwp_statespace(fwp_statespace& fwp_states, run_params run_parameters);

#endif // MAKE_WP_STATES_H
#ifndef MAKE_WP_STATES_H
#define MAKE_WP_STATES_H

#include <iostream>
#include <math.h>

#include "constants.h"
#include "error_management.h"
#include "store_functionality.h"
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

void make_p_bin_grid(int Np_WP, double* p_WP_array, run_params run_parameters);
void make_q_bin_grid(int Nq_WP, double* q_WP_array, run_params run_parameters);

void make_p_bin_quadrature_grids(int Np_WP, double* p_WP_array,
								 int Np_per_WP, double* p_array, double* wp_array);

void make_q_bin_quadrature_grids(int Nq_WP, double* q_WP_array,
								 int Nq_per_WP, double* q_array, double* wq_array);

#endif // MAKE_WP_STATES_H
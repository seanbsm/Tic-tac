#ifndef MAKE_WP_STATES_H
#define MAKE_WP_STATES_H

#include <iostream>

#include "error_management.h"

double q_normalisation(double p0, double p1);
double p_normalisation(double q0, double q1);

double p_weight_function(double p);
double q_weight_function(double q);

void make_p_bin_quadrature_grids();
void make_q_bin_quadrature_grids();

#endif // MAKE_WP_STATES_H
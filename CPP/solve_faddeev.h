#ifndef SOLVE_FADDEEV_H
#define SOLVE_FADDEEV_H

#include <iostream>
#include <math.h>
#include <complex>
#include <ctime>
#include <chrono>
#include "mkl.h"

#include "constants.h"
#include "type_defs.h"
#include "error_management.h"

void iterate_faddeev_equations(cdouble* U_array,
                               cdouble* G_array,
                               double* P123_array,
                               double* C_WP_unco_array,
					           double* C_WP_coup_array,
					           double* V_WP_unco_array,
                               double* V_WP_coup_array,
                               int Nq_WP, double* q_WP_array,
					           int Np_WP, double* p_WP_array,
					           int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array);

/* Solves the Faddeev equations
 * U = P*V + P*V*G*U
 * on the form L*U = R, where L and R are the left-
 * and right-handed sides of the equations, given by 
 * L = 1 - P*V*G
 * R = P*V
 * Since G is expressed in an SWP basis, we also must include the basis-transormation matrices C */
void direct_solve_faddeev_equations(cdouble* U_array,
                                    cdouble* G_array,
                                    double* P123_array,
                                    double* C_WP_unco_array,
					                double* C_WP_coup_array,
					                double* V_WP_unco_array,
                                    double* V_WP_coup_array,
                                    int Nq_WP,
					                int Np_WP,
					                int Nalpha,
                                    int* L_2N_array,
                                    int* S_2N_array,
                                    int* J_2N_array,
                                    int* T_2N_array,
                                    int* L_1N_array, 
                                    int* two_J_1N_array,
                                    int* two_J_3N_array,
                                    int* two_T_3N_array);

#endif // SOLVE_FADDEEV_H
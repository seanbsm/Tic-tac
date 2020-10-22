#ifndef LIPPMANN_SCHWINGER_SOLVER_H
#define LIPPMANN_SCHWINGER_SOLVER_H

#include <iostream>
#include <complex>
#include "mkl.h"

#include <string>
#include <vector>

#include "error_management.h"
#include "constants.h"
#include "type_defs.h"
#include "Interactions/potential_model.h"

void cdot_VS_colMaj(std::complex<float>*  A, std::complex<float>*  B, int N);
void cdot_VS_colMaj(std::complex<double>* A, std::complex<double>* B, int N);
void solve_MM_lin_eq(std::complex<float>*  A, std::complex<float>*  B, int N);
void solve_MM_lin_eq(std::complex<double>* A, std::complex<double>* B, int N);

//double extract_potential_element_from_array(int L, int Lp, int J, int S, bool coupled, double* V_array);

void make_denominator_array(cfloatType* D_array, int Nk, double* k_array, double* wk_array, double q, double M);

void make_wave_matrix(cfloatType* F_array, cfloatType* D_array, int Nk1, bool coupled);

double calculate_t_element(double* V_prestored_array,
						   int L, int L_p, int S, int J, int T,
						   double E, double M,
						   int Nk, double* k_array, double* wk_array,
						   int idx_row, int idx_col,
						   potential_model* pot_ptr_nn,
                           potential_model* pot_ptr_np);

#endif // LIPPMANN_SCHWINGER_SOLVER_H
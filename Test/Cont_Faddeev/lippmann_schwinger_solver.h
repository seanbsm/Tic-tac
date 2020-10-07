#ifndef LIPPMANN_SCHWINGER_SOLVER_H
#define LIPPMANN_SCHWINGER_SOLVER_H

#include <iostream>
#include <string>
#include <vector>
#include <complex>
#include "mkl.h"

#include "type_defs.h"

void cdot_VS_colMaj(std::complex<float> *A, std::complex<float> *B, int N);
void cdot_VS_colMaj(std::complex<double> *A, std::complex<double> *B, int N);

void make_denominator_array(cfloatType *D);

void make_wave_matrix(cfloatType *F, cfloatType *D, bool isCoupled);

void solve_for_T_element();

#endif // LIPPMANN_SCHWINGER_SOLVER_H
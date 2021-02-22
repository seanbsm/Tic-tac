#ifndef STORE_FUNCTIONALITY_H
#define STORE_FUNCTIONALITY_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <cstring> //std::strcpy

#include <math.h>
#include <complex>
#include "mkl.h"

#include <fstream>
#include <time.h>

#include <omp.h>

/* gsl-libraries */
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf.h>

/* hdf5-storage libaries */
#include "hdf5/serial/hdf5.h"
#include "hdf5/serial/hdf5_hl.h"

#include "error_management.h"
#include "General_functions/matrix_routines.h"

void read_sparse_matrix_elements_P123_h5 (double** P123_sparse_ptr_val_array,
                                          int**    P123_sparse_ptr_idx_array,
                                          int*     P123_sparse_ptr_dim_array,
                                          int  Np_3N, double *p_3N,
                                          int  Nq_3N, double *q_3N,
                                          int  N_chn_3N,
                                          int* chn_3N_idx_array,
                                          int  Jj_dim,
                                          int* L12_Jj,
                                          int* S12_Jj,
                                          int* J12_Jj,
                                          int* T12_Jj,
                                          int* l3_Jj,
                                          int* two_j3_Jj,
                                          int* two_J_Jj,
                                          int* two_T_Jj,
                                          int* P_3N_array,
                                          std::string filename_in);

void store_sparse_matrix_elements_P123_h5 (double** P123_sparse_ptr_val_array,
                                           int**    P123_sparse_ptr_idx_array,
                                           int*     P123_sparse_ptr_dim_array,
                                           int  Np_3N, double *p_3N, int Nq_3N, double *q_3N,
                                           int  N_chn_3N,
                                           int* chn_3N_idx_array,
                                           int  Jj_dim,
                                           int *L12_Jj,
                                           int *S12_Jj,
                                           int *J12_Jj,
                                           int *T12_Jj,
                                           int *l3_Jj,
                                           int *two_j3_Jj,
                                           int *two_J_Jj,
                                           int *two_T_Jj,
                                           int *P_3N_array,
                                           std::string filename_in);

void store_matrix_elements_P123_h5 (double *P123_store,
                                    int Np_3N, double *p_3N, int Nq_3N, double *q_3N,
                                    int Jj_dim,
                                    int *L12_Jj,
                                    int *S12_Jj,
                                    int *J12_Jj,
                                    int *T12_Jj,
                                    int *l3_Jj,
                                    int *two_j3_Jj,
                                    int *two_J_Jj,
                                    int *two_T_Jj,
                                    std::string filename_in);

#endif // STORE_FUNCTIONALITY_H
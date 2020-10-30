#ifndef PERMUTATION_OPERATORS_H
#define PERMUTATION_OPERATORS_H

#include <string>
#include <cstring> //std::strcpy
#include <iostream>
#include <fstream>
#include <complex>
#include <algorithm>
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream
#include <vector>

/* Modules used in benchmark code */
#include "mkl.h"
#include <math.h>
#include "auxiliary.h"

/* GSL functionality */
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf.h>

/* .h5-file functionality */
#include "hdf5/serial/hdf5.h"
#include "hdf5/serial/hdf5_hl.h"

/* Project headers */
#include "error_management.h"

typedef struct pw_table
{
    int alpha_pwtable;
    int L_pwtable;
    int S_pwtable;
    int J_pwtable;
    int T_pwtable;
    int l_pwtable;
    int twoj_pwtable;
    int twoJtotal_pwtable;
    int PARtotal_pwtable;
    int twoTtotal_pwtable;
} pw_table;

typedef struct alpha_table
{
    int alpha_alphaprime_index;
    int alpha_index;
    int alphaprime_index;
} alpha_table;


typedef struct pmesh_table
{
    int index_ptable;
    double mesh_ptable;
    double weight_ptable;
} pmesh_table;


typedef struct qmesh_table
{
    int index_qtable;
    double mesh_qtable;
    double weight_qtable;
} qmesh_table;

void check_h5_read_call(herr_t ret);
void check_h5_read_table_call(herr_t ret);
void check_h5_close_call(herr_t ret);

void read_P123_h5_data_file(std::string file_path, double* P123, int Nq_3N, double *q_3N, int Np_3N, double *p_3N);

void get_h5_P123_dimensions(std::string file_path, int& Nalpha, int& Np, int& Nq);

void calculate_permutation_operator(double* P123_array,
                                    int Nq, double* q_array, double* wq_array,
                                    int Np, double* p_array, double* wp_array,
                                    int Nx, double* x_array, double* wx_array,
                                    int Nalpha,
                                    int* L_2N_array,
                                    int* S_2N_array,
                                    int* J_2N_array,
                                    int* T_2N_array,
                                    int* l_3N_array,
                                    int* two_j_3N_array,
                                    int two_J_3N, int two_T_3N, int parity_3N);

void calculate_antisymmetrization_operator(std::string file_path, 
                                           int &Np, int &Nq, int& Nalpha,
                                           double** A123,
                                           double* q_array,
                                           double* p_array);

#endif // PERMUTATION_OPERATORS_H
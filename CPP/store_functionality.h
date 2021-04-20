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

void store_sparse_permutation_matrix_for_3N_channel_h5(double* P123_sparse_val_array,
													   int*    P123_sparse_row_array,
													   int*    P123_sparse_col_array,
													   size_t  P123_sparse_dim,
													   int     Np_WP, double* p_WP_array,
													   int     Nq_WP, double* q_WP_array,
													   int     Nalpha,
													   int*    L_2N_array,
													   int*    S_2N_array,
													   int*    J_2N_array,
													   int*    T_2N_array,
													   int*    L_1N_array, 
													   int*    two_J_1N_array,
													   int     two_J_3N,
													   int     two_T_3N,
													   int     P_3N,
													   std::string filename_in,
													   bool print_content);
													   
void read_sparse_permutation_matrix_for_3N_channel_h5(double** P123_sparse_val_array,
													   int**   P123_sparse_row_array,
													   int**   P123_sparse_col_array,
													   size_t& P123_sparse_dim,
													   int     Np_WP, double* p_WP_array,
													   int     Nq_WP, double* q_WP_array,
													   int     Nalpha,
													   int*    L_2N_array,
													   int*    S_2N_array,
													   int*    J_2N_array,
													   int*    T_2N_array,
													   int*    L_1N_array, 
													   int*    two_J_1N_array,
													   int     two_J_3N,
													   int     two_T_3N,
													   int     P_3N,
													   std::string filename_in,
													   bool print_content);

void write_integer_to_h5(int integer, char* int_name, hid_t file_id);
void read_integer_from_h5(int& integer, char* int_name, char* filename);

void write_ULL_integer_to_h5(unsigned long long int integer, char* int_name, hid_t file_id);
void read_ULL_integer_from_h5(unsigned long long& integer, char* int_name, char* filename);

void write_WP_boundaries_to_h5(double* WP_boundaries, int N_WP, char* mesh_name, hid_t file_id);
void read_WP_boundaries_from_h5(double* WP_boundaries, int N_WP, char* mesh_name, char* filename);

void write_PW_statespace_to_h5(int   Nalpha,
							   int*  L_2N_array,
							   int*  S_2N_array,
							   int*  J_2N_array,
							   int*  T_2N_array,
							   int*  L_1N_array, 
							   int*  two_J_1N_array,
							   int   two_J_3N,
							   int   two_T_3N,
							   int   P_3N,
							   hid_t file_id);
void read_PW_statespace_to_h5(int  Nalpha,
							  int* L_2N_array,
							  int* S_2N_array,
							  int* J_2N_array,
							  int* T_2N_array,
							  int* L_1N_array, 
							  int* two_J_1N_array,
							  int* two_J_3N_array,
							  int* two_T_3N_array,
							  int* P_3N_array,
							  char* filename);

void write_sparse_permutation_matrix_h5(double* P123_sparse_val_array,
										int*    P123_sparse_row_array,
										int*    P123_sparse_col_array,
										size_t  P123_sparse_dim,
										hid_t   file_id);
void read_sparse_permutation_matrix_h5(double* P123_sparse_val_array,
									   int*    P123_sparse_row_array,
									   int*    P123_sparse_col_array,
									   size_t  P123_sparse_dim,
									   char*   filename);

#endif // STORE_FUNCTIONALITY_H
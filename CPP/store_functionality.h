#ifndef STORE_FUNCTIONALITY_H
#define STORE_FUNCTIONALITY_H

#include <iostream>
#include <iomanip>
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

#include "type_defs.h"
#include "error_management.h"
#include "store_functionality.h"
#include "General_functions/matrix_routines.h"

void store_sep_complex_matrix(double* re_matrix_array,
							  double* im_matrix_array,
                         	  size_t num_rows,
						 	  size_t num_cols,
						 	  size_t steplength,
						 	  std::string file_path,
						 	  bool rewrite_file,
						 	  std::string array_seperator_text);

void store_complex_matrix(cdouble* matrix_array,
                          size_t num_rows,
						  size_t num_cols,
						  size_t steplength,
						  std::string file_path,
						  bool rewrite_file,
						  std::string array_seperator_text);

void store_run_parameters(run_params run_parameters);

void store_q_WP_kinematics_txt(size_t Nq_WP,
							   double* q_WP_array,
							   double* Eq_WP_array,
							   double* Tlab_WP_array,
							   double* q_WP_midpoint_array,
							   double* Eq_WP_midpoint_array,
							   double* Tlab_WP_midpoint_array,
						       std::string filename);
void store_q_WP_boundaries_csv(size_t Nq_WP, double* q_WP_array,
						   	   std::string filename);
void store_U_matrix_elements_txt(std::complex<double>*    U_array,
								 std::string 			  potential_model,
								 int					  two_J,
								 int					  P_3N,
								 int					  Np_WP,
								 int					  Nq_WP,
								 double					  E_bound,
								 double*				  T_lab_array,
								 double*				  E_com_array,
							     int* q_com_idx_array,	  size_t num_q_com,
					  		     int* deuteron_idx_array, size_t num_deuteron_states,
							     int* L_1N_array, 
							     int* two_J_1N_array,
							     std::string filename);
void store_U_matrix_elements_csv(std::complex<double>* U_array,
							 	 int* q_com_idx_array,	  size_t num_q_com,
					  		 	 int* deuteron_idx_array, size_t num_deuteron_states,
							 	 int* L_1N_array, 
							 	 int* two_J_1N_array,
							 	 std::string filename);

void read_WP_boundaries_from_txt(double* WP_boundaries,
								 int N_WP,
								 std::string filename);

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
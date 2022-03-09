#ifndef DISK_IO_FUNCTIONALITY_H
#define DISK_IO_FUNCTIONALITY_H

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
#include "set_run_parameters.h"
#include "error_management.h"
#include "General_functions/matrix_routines.h"

void read_parameter_sample_list(std::string filename,
								std::vector<double>& parameter_vector,
								int& num_param,
								int& num_param_sets);

void read_input_energies(double*& energy_array,
						 int&	  num_energies,
						 std::string file_path);

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
void store_q_WP_boundaries_csv(fwp_statespace fwp_states,
						   	   std::string filename);
void store_U_matrix_elements_txt(std::complex<double>*    U_array,
								 solution_configuration solve_config,
								 channel_os_indexing chn_os_indexing,
								 run_params run_parameters,
								 swp_statespace swp_states,
								 pw_3N_statespace pw_states,
							     std::string filename);
void store_U_matrix_elements_csv(std::complex<double>* U_array,
							 	 int* q_com_idx_array,	  size_t num_q_com,
					  		 	 int* deuteron_idx_array, size_t num_deuteron_states,
							     pw_3N_statespace pw_states,
							 	 std::string filename);

void read_WP_boundaries_from_txt(double* WP_boundaries,
								 int N_WP,
								 std::string filename);

void store_sparse_matrix_h5(double* sparse_val_array,
							int*    sparse_row_array,
							int*    sparse_col_array,
							size_t  sparse_dim,
							size_t  dense_dim,
							std::string filename_in,
							bool print_content);
void read_sparse_matrix_h5(double** sparse_val_array,
						   int**    sparse_row_array,
						   int**    sparse_col_array,
						   size_t&  sparse_dim,
						   size_t&  dense_dim,
						   std::string filename_in,
						   bool    print_content);

void store_sparse_permutation_matrix_for_3N_channel_h5(double* P123_sparse_val_array,
													   int*    P123_sparse_row_array,
													   int*    P123_sparse_col_array,
													   size_t  P123_sparse_dim,
													   int     Np_WP, double* p_WP_array,
													   int     Nq_WP, double* q_WP_array,
													   pw_3N_statespace pw_states,
													   std::string filename_in,
													   bool print_content);
													   
void read_sparse_permutation_matrix_for_3N_channel_h5(double** P123_sparse_val_array,
													   int**   P123_sparse_row_array,
													   int**   P123_sparse_col_array,
													   size_t& P123_sparse_dim,
													   int     Np_WP, double* p_WP_array,
													   int     Nq_WP, double* q_WP_array,
													   pw_3N_statespace pw_states,
													   std::string filename_in,
													   bool print_content);

void write_integer_to_h5(int integer, const char* int_name, hid_t file_id);
void read_integer_from_h5(int& integer, const char* int_name, const char* filename);

void write_ULL_integer_to_h5(unsigned long long int integer, const char* int_name, hid_t file_id);
void read_ULL_integer_from_h5(unsigned long long int& integer, const char* int_name, const char* filename);

void write_WP_boundaries_to_h5(double* WP_boundaries, int N_WP, const char* mesh_name, hid_t file_id);
void read_WP_boundaries_from_h5(double* WP_boundaries, int N_WP, const char* mesh_name, const char* filename);

void write_PW_statespace_to_h5(int   Nalpha,
							   int*  L_2N_array,
							   int*  S_2N_array,
							   int*  J_2N_array,
							   int*  T_2N_array,
							   int*  L_1N_array, 
							   int*  two_J_1N_array,
							   int*  two_T_3N_array,
							   int   two_J_3N,
							   int   P_3N,
							   hid_t file_id);
void read_PW_statespace_to_h5(int  Nalpha,
							  int* L_2N_array,
							  int* S_2N_array,
							  int* J_2N_array,
							  int* T_2N_array,
							  int* L_1N_array, 
							  int* two_J_1N_array,
							  int* two_T_3N_array,
							  int* two_J_3N_array,
							  int* P_3N_array,
							  const char* filename);

void write_sparse_permutation_matrix_h5(double* P123_sparse_val_array,
										int*    P123_sparse_row_array,
										int*    P123_sparse_col_array,
										size_t  P123_sparse_dim,
										hid_t   file_id);
void read_sparse_permutation_matrix_h5(double* P123_sparse_val_array,
									   int*    P123_sparse_row_array,
									   int*    P123_sparse_col_array,
									   size_t  P123_sparse_dim,
									   const char*   filename);

#endif // DISK_IO_FUNCTIONALITY_H
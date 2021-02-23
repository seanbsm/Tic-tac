
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>

/* Used for printf */
#include <cstdio>
#include <cstring>

/* Time-keeping modules */
#include <chrono>
#include <ctime>

#include "make_pw_symm_states.h"
#include "make_permutation_matrix.h"
#include "General_functions/gauss_legendre.h"
#include "Interactions/potential_model.h"
#include "make_potential_matrix.h"
#include "make_swp_states.h"
#include "make_resolvent.h"
#include "solve_faddeev.h"
#include "basis_transformations.h"
#include "store_functionality.h"

using namespace std;

int main(int argc, char* argv[]){

	auto program_start = chrono::system_clock::now();
	
	/* ------------------- Start main body of code here ------------------- */
	/* Start of code segment for parameters, variables and arrays declaration */

	/* Current scattering energy */
	double E = 10;

	/* Setting to store calculated P123 matrix in WP basis to h5-file */
	bool calculate_and_store_P123 = true;

	/* PWE truncation */
	/* Maximum (max) values for J_2N and J_3N (minimum is set to 0 and 1, respectively)*/
    int J_2N_max 	 = 0;//5;
	int two_J_3N_max = 1;//25;
	if ( two_J_3N_max%2==0 ||  two_J_3N_max<=0 ){
		raise_error("Cannot have even two_J_3N_max");
	}

	/* Wave-packet 3N momenta */
	int Np_WP	   	 = 40;
	int Nq_WP	   	 = 40;
	double* p_WP_array  = NULL;
	double* q_WP_array  = NULL;

	/* Quadrature 3N momenta per WP cell */
	int Np_per_WP	 = 8;
	int Nq_per_WP	 = 8;
	double* p_array  = NULL;
	double* q_array  = NULL;
	double* wp_array = NULL;
	double* wq_array = NULL;

	/* Permutation matrix - read from file */
	string  P123_file_path = "../../Data/3N_permutation_operator/P123_files/P123_medium.h5";
	double* P123_array    = NULL;
	double* P123_p_array  = NULL;
	double* P123_q_array  = NULL;
	double* P123_wp_array = NULL;
	double* P123_wq_array = NULL;

	/* Potential matrices */
	double* V_WP_unco_array = NULL;
    double* V_WP_coup_array = NULL;

	/* Quantum numbers of partial-wave expansion in state_3N_array.
	 * All non-specified quantum numbers are either given by nature,
	 * deduced from the values below, or disappear in summations later */
	int  Nalpha   		= 0;		// Number of partial waves, set in dynamical state-space construction (by construct_symmetric_pw_states)
    int* L_2N_array     = NULL;     // pair-nucleon state angular momentum
    int* S_2N_array     = NULL;     // pair-nucleon state total spin
    int* J_2N_array     = NULL;     // pair-nucleon state total angular momentum
    int* T_2N_array     = NULL;     // pair-nucleon state total isospin
    int* L_1N_array     = NULL;     // orbital-nucleon state angular momentum
    int* two_J_1N_array = NULL; 	// orbital-nucleon state total angular momentum x2
	int* two_J_3N_array = NULL; 	// three-nucleon state total angular momentum x2
	int* two_T_3N_array = NULL; 	// three-nucleon state total isospin x2
	int* P_3N_array 	= NULL; 	// three-nucleon state parity

	/* Three-nucleon channel indexing (a 3N channel is defined by a set {two_J_3N, two_T_3N, P_3N}) */
	int  N_chn_3N 		  = 0;
	int* chn_3N_idx_array = NULL;


	/* Potential model class pointers */
	potential_model* pot_ptr_np = NULL;
	potential_model* pot_ptr_nn = NULL;

	/* printf buffer used to add text to current terminal printout */
	//char buffer[1024];
	//setvbuf(stdout, buffer, _IOFBF, 1024);

	/* End of code segment for variables and arrays declaration */
	/* Start of code segment for state space construction */

	printf("Constructing 3N partial-wave basis ... \n");
	construct_symmetric_pw_states(J_2N_max,
								  two_J_3N_max,
                                  N_chn_3N,
                                  &chn_3N_idx_array,
								  Nalpha,
								  &L_2N_array,
								  &S_2N_array,
								  &J_2N_array,
								  &T_2N_array,
								  &L_1N_array,
								  &two_J_1N_array,
								  &two_J_3N_array,
								  &two_T_3N_array,
								  &P_3N_array);
	printf(" - Done \n");

	printf("Constructing wave-packet (WP) p-momentum bin boundaries ... \n");
	p_WP_array = new double [Np_WP+1];
	make_p_bin_grid(Np_WP, p_WP_array);
	printf(" - Done \n");
	printf("Constructing wave-packet (WP) q-momentum bin boundaries ... \n");
	q_WP_array = new double [Nq_WP+1];
	make_q_bin_grid(Nq_WP, q_WP_array);
	printf(" - Done \n");

	printf("Constructing p quadrature mesh per WP, for all WPs ... \n");
	p_array  = new double [Np_per_WP*Np_WP];
	wp_array = new double [Np_per_WP*Np_WP];
	make_p_bin_quadrature_grids(Np_WP, p_WP_array,
                                Np_per_WP, p_array, wp_array);
	printf(" - Done \n");
	printf("Constructing q quadrature mesh per WP, for all WPs ... \n");
	q_array  = new double [Nq_per_WP*Nq_WP];
	wq_array = new double [Nq_per_WP*Nq_WP];
	make_q_bin_quadrature_grids(Nq_WP, q_WP_array,
                                Nq_per_WP, q_array, wq_array);
	printf(" - Done \n");

	/* End of code segment for state space construction */
	/* Start of code segment for permutation matrix construction */
		
	int P123_dim = Nalpha*Np_WP*Nq_WP;
	P123_array = NULL;//new double [P123_dim * P123_dim];

	double** P123_sparse_ptr_val_array = new double* [N_chn_3N];
	int** 	 P123_sparse_ptr_row_array = new int* 	 [N_chn_3N];
	int** 	 P123_sparse_ptr_col_array = new int* 	 [N_chn_3N];
	int*	 P123_sparse_ptr_dim_array = new int 	 [N_chn_3N];
		
	if (calculate_and_store_P123){

		int Nx 			 = 20;
		double* x_array  = new double [Nx];
    	double* wx_array = new double [Nx];
		gauss(x_array, wx_array, Nx);

		printf("Calculating P123 ... \n");
		auto timestamp_P123_calc_start = chrono::system_clock::now();
		calculate_permutation_matrices_for_all_3N_channels(P123_sparse_ptr_val_array,
    	                                            	   P123_sparse_ptr_row_array,
    	                                            	   P123_sparse_ptr_col_array,
    	                                            	   P123_sparse_ptr_dim_array,
    	                                            	   Nq_WP*Nq_per_WP, q_array, wq_array, Np_per_WP, Np_WP, p_WP_array,
    						                           	   Np_WP*Np_per_WP, p_array, wp_array, Nq_per_WP, Nq_WP, q_WP_array,
    						                           	   Nx, x_array, wx_array,
														   N_chn_3N,
    	                            					   chn_3N_idx_array,
    						                           	   Nalpha,
    						                           	   L_2N_array,
    						                           	   S_2N_array,
    						                           	   J_2N_array,
    						                           	   T_2N_array,
    						                           	   L_1N_array,
    						                           	   two_J_1N_array,
    						                         	   two_J_3N_array,
    						                         	   two_T_3N_array);
		auto timestamp_P123_calc_end = chrono::system_clock::now();
		chrono::duration<double> time_P123_calc = timestamp_P123_calc_end - timestamp_P123_calc_start;
		printf(" - Done. Time used: %.6f\n", time_P123_calc.count());

		printf("Storing P123 to h5 ... \n");
		auto timestamp_P123_store_start = chrono::system_clock::now();
		store_sparse_matrix_elements_P123_h5 (P123_sparse_ptr_val_array,
    	                            		  P123_sparse_ptr_row_array,
    	                                      P123_sparse_ptr_col_array,
    	                            		  P123_sparse_ptr_dim_array,
    	                            		  Np_WP, p_WP_array, Nq_WP, q_WP_array,
											  N_chn_3N,
    	                            		  chn_3N_idx_array,
    	                            		  Nalpha,
    	                            		  L_2N_array,
    	                          			  S_2N_array,
    	                          			  J_2N_array,
    	                          			  T_2N_array,
    	                          			  L_1N_array,
    	                          			  two_J_1N_array,
    	                        			  two_J_3N_array,
    	                        			  two_T_3N_array,
											  P_3N_array,
    	                            		  "P123_sparse.h5");
		auto timestamp_P123_store_end = chrono::system_clock::now();
		chrono::duration<double> time_P123_store = timestamp_P123_store_end - timestamp_P123_store_start;
		printf(" - Done. Time used: %.6f\n", time_P123_store.count());
	}
	else{
		printf("Reading P123 from h5 ... \n");

		double** P123_sparse_ptr_val_array_t = new double* [N_chn_3N];
		int** 	 P123_sparse_ptr_row_array_t = new int* 	 [N_chn_3N];
		int** 	 P123_sparse_ptr_col_array_t = new int* 	 [N_chn_3N];
		int*	 P123_sparse_ptr_dim_array_t = new int 	 [N_chn_3N];

		auto timestamp_P123_read_start = chrono::system_clock::now();
		read_sparse_matrix_elements_P123_h5 (P123_sparse_ptr_val_array_t,
    	                            		 P123_sparse_ptr_row_array_t,
    	                            		 P123_sparse_ptr_col_array_t,
    	                            		 P123_sparse_ptr_dim_array_t,
    	                            		 Np_WP+1, p_WP_array, Nq_WP+1, q_WP_array,
											 N_chn_3N,
    	                            		 chn_3N_idx_array,
    	                            		 Nalpha,
    	                            		 L_2N_array,
    	                          			 S_2N_array,
    	                          			 J_2N_array,
    	                          			 T_2N_array,
    	                          			 L_1N_array,
    	                          			 two_J_1N_array,
    	                        			 two_J_3N_array,
    	                        			 two_T_3N_array,
											 P_3N_array,
    	                            		 "P123_sparse.h5");
		auto timestamp_P123_read_end = chrono::system_clock::now();
		chrono::duration<double> time_P123_read = timestamp_P123_read_end - timestamp_P123_read_start;
		printf(" - Done. Time used: %.6f\n", time_P123_read.count());
	}

	/* End of code segment for permutation matrix construction */
	/* Start of code segment for potential matrix construction */

	//pot_ptr_np = potential_model::fetch_potential_ptr("LO_internal", "np");
	//pot_ptr_nn = potential_model::fetch_potential_ptr("LO_internal", "nn");
	//pot_ptr_np = potential_model::fetch_potential_ptr("N2LOopt", "np");
	//pot_ptr_nn = potential_model::fetch_potential_ptr("N2LOopt", "nn");
	pot_ptr_np = potential_model::fetch_potential_ptr("Idaho_N3LO", "np");
	pot_ptr_nn = potential_model::fetch_potential_ptr("Idaho_N3LO", "nn");

	printf("Constructing 2N-potential matrices in WP basis ... \n");
	int V_unco_array_size = Np_WP*Np_WP   * 2*(J_2N_max+1);
    int V_coup_array_size = Np_WP*Np_WP*4 *    J_2N_max;
    V_WP_unco_array = new double [V_unco_array_size];
    V_WP_coup_array = new double [V_coup_array_size];
	calculate_potential_matrices_array_in_WP_basis(V_WP_unco_array,
                                                   V_WP_coup_array,
												   true,
                                                   Np_WP, p_WP_array,
                                                   Np_per_WP, p_array, wp_array,
                                                   Nalpha, L_2N_array, S_2N_array, J_2N_array, T_2N_array,
                                                   pot_ptr_nn,
                                                   pot_ptr_np);
	printf(" - Done \n");

	/* End of code segment for potential matrix construction */
	/* Start of code segment for scattering wave-packets construction */

	double* p_SWP_unco_array = new double [  (Np_WP+1) * 2*(J_2N_max+1)];
	double* p_SWP_coup_array = new double [2*(Np_WP+1) *    J_2N_max];

	double* C_WP_unco_array = new double [V_unco_array_size];
	double* C_WP_coup_array = new double [V_coup_array_size];
	
	printf("Constructing 2N SWPs ... \n");
	make_swp_states(p_SWP_unco_array,
					p_SWP_coup_array,
					C_WP_unco_array,
					C_WP_coup_array,
					V_WP_unco_array,
                    V_WP_coup_array,
					Np_WP, p_WP_array,
					Nalpha, L_2N_array, S_2N_array, J_2N_array, T_2N_array,
					J_2N_max);
	printf(" - Done \n");

	/* End of code segment for scattering wave-packets construction */
	/* Start of code segment for resolvent matrix (diagonal array) construction */

	/* Resolvent array */
	cdouble* G_array = new cdouble [Nalpha * Nq_WP * Np_WP];

	printf("Constructing 3N resolvents ... \n");
	calculate_resolvent_array_in_SWP_basis(G_array,
                                           E,
                                           Np_WP,
                                           p_SWP_unco_array,
					                       p_SWP_coup_array,
					                       Nq_WP, q_WP_array,
					                       Nalpha, L_2N_array, S_2N_array, J_2N_array, T_2N_array);
	printf(" - Done \n");

	/* End of code segment for resolvent matrix (diagonal array) construction */
	/* Start of code segment for iterations of elastic Faddeev equations */

	cdouble* U_array = new cdouble [Nalpha*Np_WP*Nq_WP * Nalpha*Np_WP*Nq_WP];

	printf("Solving Faddeev equations ... \n");
	direct_solve_faddeev_equations(U_array,
                                   G_array,
                                   P123_array,
                                   C_WP_unco_array,
					               C_WP_coup_array,
					               V_WP_unco_array,
                                   V_WP_coup_array,
                                   Nq_WP,
					               Np_WP,
					               Nalpha,
                                   L_2N_array,
                                   S_2N_array,
                                   J_2N_array,
                                   T_2N_array,
                                   L_1N_array, 
                                   two_J_1N_array,
								   two_J_3N_array,
								   two_T_3N_array);
	printf(" - Done \n");

	/* End of code segment for iterations of elastic Faddeev equations */
	/* Start of code segment for calculating scattering observables */


	/* End of code segment for calculating scattering observables */


	/* -------------------- End main body of code here -------------------- */

	auto program_end= chrono::system_clock::now();

	chrono::duration<double> total_time = program_end - program_start;
	printf("Total run-time: %.6f\n", total_time.count());

	printf("END OF RUN");
	
    return 0;
}


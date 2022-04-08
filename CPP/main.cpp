
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

#include "set_run_parameters.h"
#include "make_pw_symm_states.h"
#include "make_permutation_matrix.h"
#include "General_functions/gauss_legendre.h"
#include "General_functions/kinetic_conversion.h"
#include "Interactions/potential_model.h"
#include "make_potential_matrix.h"
#include "make_wp_states.h"
#include "make_swp_states.h"
#include "make_resolvent.h"
#include "solve_faddeev.h"
#include "disk_io_routines.h"
#include "run_organizer.h"

/*

Faddev eqs.
U = C^T PVC + C^T PVC G U

PROGRAM STRUCTURE:
Define PW state space:
 - input:  J_2N_max, J_3N_max, T_3N_max
 - output: {alpha=(J,L,S,T,l,j,J3,T3,P3)} (9 quantum numbers per state)

Define free WP space:
 - input: grid type (chebyshev), Np_WP, Nq_WP, sparness degree (t), scaling, momentum/energy WP
 - output: {q} and {p} WP boundaries (cells)

Calculate 2N interations:
 - input: {alpha}->{n}, {p} (basically a 2N-basis loop for J_2N_max)
 - output: V_2N_uncoupled_array and V_2N_coupled_array

Construct SWP space from Hamitonian diagonalisation of 2N potentials:
 - input: {alpha}->{n}, {p}, V_2N_... (basically a 2N-basis loop for J_2N_max)
 - output: C_2N_uncoupled_array, C_2N_coupled_array, {p}_SWP-boundaries for SWPs (which are {n}-dependent)

For-loop over 3N channels chn"="(J3,T3,P3) ({alpha}_chn stored congruently in memory):
 - Calculate/read P123:
  - input: {alpha}_chn, {p}, {q}
  - output: sparse P123 matrix P123_array in COO format (COO- to CSR-format converter available)
  - comment: stored to disk if boolean "calculate" is chosen

 - Calculate 3N resolvent from SWP {p}-boundaries:
   - input:  {alpha}_chn, {p}_SWP, {q}, {E} (on-shell energies)
   - output: G_array of length len{alpha} x len{p} x len{q} x len{E}
 
 - Call function run_parameters.solve_faddeev_equation (Nalpha_chn = number of alpha in current 3N channel)
   - input: Nalpha_chn, Np_WP, Nq_WP, V_2N_..., C_2N_..., G_array, P123_array
   - method (Pade): U = A + KU; A = C^T PVC; K=C^T PVC G = AG
			          = A + KA + KKA + KKKA + KKKKA + .... (Neumann sum)
			          = A + AGA + AGAGA + AGAGAGA + ...
			          = sum_n=0^N A(GA)^n
			          = sum_n=0^N (AG)^n A
	- method (DSS): Rewrite Faddeev: (1-AG)U = A_c, where A_c is a column of interest (on-shell)
	- bottleneck: A=C^T PVC column-calculation
	- output: Elastic U-matrix elements (several {n}=deutron channels within one "chn")
	
A^n = (C^T PVC)(C^T PVC)(C^T PVC)(C^T PVC)...

A   =      C^t PV C
A^2 =  A   C^t PV C = C^t PV C C^t PV C          = C^t PVPV C
A^3 =  A^2 C^t PV C = C^t PV C C^t PV C C^t PV C = C^t PVPVPV C

*/

using namespace std;

int main(int argc, char* argv[]){

	auto program_start = chrono::system_clock::now();
	
	// Handy for converting Tlab to Ecm
	//std::cout << com_q_momentum_to_com_energy(lab_energy_to_com_momentum(13, -Ed_measured)) << std::endl;
	//return 0;

	/* ------------------- Start main body of code here ------------------- */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* Start of code segment for parameters, variables and arrays declaration */

	/* Interpret command-line input and save in run_parameters */
	run_params run_parameters;
	set_run_parameters(argc, argv, run_parameters);

	/* Wave-packet 3N momenta */
	int Np_WP = run_parameters.Np_WP;
	int Nq_WP = run_parameters.Nq_WP;

	/* PWE truncation */
	/* Maximum (max) values for J_2N and J_3N (minimum is set to 0 and 1, respectively)*/
	int J_2N_max 	 = run_parameters.J_2N_max;

	/* Quadrature 3N momenta per WP cell */
	fwp_statespace fwp_states;

	/* Quantum numbers of partial-wave expansion in state_3N_array.
	 * All non-specified quantum numbers are either given by nature,
	 * deduced from the values below, or disappear in summations later */
	pw_3N_statespace pw_states;

	/* Potential model class pointers */
	potential_model* pot_ptr  = potential_model::fetch_potential_ptr(run_parameters);

	/* End of code segment for variables and arrays declaration */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* Start of code segment for reading input model parameters (if enabled) */

	std::vector<double> parameter_vector;
	int num_params_to_loop = run_parameters.PSI_end - run_parameters.PSI_start;
	int	num_params_in_file;
	int num_model_params;
	read_parameter_sample_list(run_parameters, parameter_vector, num_model_params, num_params_in_file);
	if (num_params_to_loop>num_params_in_file){
		raise_error("PSI-range given is larger than parameter-set input file.");
	}
	//std::cout << num_param_sets << " " << num_params << std::endl;
	//for (const auto &val : parameter_vector){
	//	std::cout << val << std::endl;
	//}

	/* End of code segment for reading input model parameters (if enabled) */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* Start of code segment for state space construction */
	
	printf("Constructing 3N partial-wave basis ... \n");
	construct_symmetric_pw_states(pw_states,
								  run_parameters);

	/* Allocate deuteron-channel index-lookup arrays */
	solution_configuration solve_config;
	solve_config.deuteron_idx_arrays = new int* [pw_states.N_chn_3N];
	solve_config.deuteron_num_array  = new int  [pw_states.N_chn_3N];
	
	make_fwp_statespace(fwp_states, run_parameters);

	/* End of code segment for state space construction */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* Start of looping over 3N-channels */


	for (int chn_3N=0; chn_3N<pw_states.N_chn_3N; chn_3N++){

		if (run_parameters.parallel_run==true && chn_3N!=run_parameters.channel_idx){
			continue;
		}

		/* ################################################################################################################### */
		/* ################################################################################################################### */
		/* ################################################################################################################### */
		/* Start of 3N-channel setup */
		/* Lower and upper limits on PW state space for channel */
		int idx_alpha_lower  = pw_states.chn_3N_idx_array[chn_3N];
		int idx_alpha_upper  = pw_states.chn_3N_idx_array[chn_3N+1];

		/* Pointers to sub-arrays of PW state space corresponding to chn_3N */
		pw_3N_statespace pw_substates;
		pw_substates.Nalpha			= idx_alpha_upper - idx_alpha_lower;
		pw_substates.J_2N_max		= pw_states.J_2N_max;
		pw_substates.L_2N_array 	= &pw_states.L_2N_array[idx_alpha_lower];
		pw_substates.S_2N_array 	= &pw_states.S_2N_array[idx_alpha_lower];
		pw_substates.J_2N_array 	= &pw_states.J_2N_array[idx_alpha_lower];
		pw_substates.T_2N_array 	= &pw_states.T_2N_array[idx_alpha_lower];
		pw_substates.L_1N_array 	= &pw_states.L_1N_array[idx_alpha_lower];
		pw_substates.two_J_1N_array = &pw_states.two_J_1N_array[idx_alpha_lower];
		pw_substates.two_T_3N_array = &pw_states.two_T_3N_array[idx_alpha_lower];
		pw_substates.two_J_3N_array = &pw_states.two_J_3N_array[idx_alpha_lower];
		pw_substates.P_3N_array		= &pw_states.P_3N_array[idx_alpha_lower];
		
		int two_J_3N = pw_substates.two_J_3N_array[0];
		int P_3N 	 = pw_substates.P_3N_array[0];

		printf("Working on 3N-channel J_3N=%.d/2, PAR=%.d (channel %.d of %.d) with %.d partial-wave states \n",
				two_J_3N, P_3N, chn_3N+1, pw_states.N_chn_3N, pw_substates.Nalpha);

		/* End of 3N-channel setup */
		/* ################################################################################################################### */
		/* ################################################################################################################### */
		/* ################################################################################################################### */
		/* Start of code segment for permutation matrix construction */
		double* P123_sparse_val_array = NULL;
		int* 	P123_sparse_row_array = NULL;
		int* 	P123_sparse_col_array = NULL;
		size_t	P123_sparse_dim		  = 0;

		fill_P123_arrays(&P123_sparse_val_array,
						&P123_sparse_row_array,
						&P123_sparse_col_array,
						P123_sparse_dim,
						run_parameters.production_run,
						fwp_states,
						pw_substates,
						run_parameters,
						run_parameters.P123_folder);

		/* End of code segment for permutation matrix construction */
		/* ################################################################################################################### */
		/* ################################################################################################################### */
		/* ################################################################################################################### */
		if (run_parameters.solve_faddeev){

			/* Convert row-major sparse format to column-major */
			printf(" - Converting P123 from COO to CSC format ... \n");
			printf("   - Converting row- to column-major format ... \n");
			size_t  dense_dim = pw_substates.Nalpha * fwp_states.Nq_WP * fwp_states.Np_WP;
			unsorted_sparse_to_coo_col_major_sorter(&P123_sparse_val_array,
													&P123_sparse_row_array,
													&P123_sparse_col_array,
													P123_sparse_dim,
													dense_dim);
			printf("     - Done \n");

			/* Convert from COO format to CSC format */
			printf("   - Converting COO to CSC format ... \n");
			size_t* P123_sparse_col_array_csc = new size_t [dense_dim+1];
			coo_to_csc_format_converter(P123_sparse_col_array,
										P123_sparse_col_array_csc,
										P123_sparse_dim,
										dense_dim);
			printf("     - Done \n");

			printf(" - Looping through input parameter sets ... \n");
			for (int idx_param_set=run_parameters.PSI_start; idx_param_set<run_parameters.PSI_end; idx_param_set++){

				printf("   - Setting model parameters for PSI %d ... \n", idx_param_set);
				if (idx_param_set==run_parameters.PSI_start){
					pot_ptr->first_parameter_sampling(true);
				}
				else{
					pot_ptr->first_parameter_sampling(false);
				}
				pot_ptr->update_parameters(&parameter_vector[idx_param_set*num_model_params]);
				printf("     - Done \n");

				/* ################################################################################################################### */
				/* ################################################################################################################### */
				/* ################################################################################################################### */
				/* Start of code segment for potential matrix construction */

				int num_2N_unco_states = 0;
				int num_2N_coup_states = 0;
				if (run_parameters.tensor_force==true){
					num_2N_unco_states = 2*(J_2N_max+1);
					num_2N_coup_states =    J_2N_max;
					if (run_parameters.isospin_breaking_1S0==true){
						num_2N_unco_states -= 1;
						num_2N_coup_states += 1;
					}
				}
				else{
					num_2N_unco_states = 4*J_2N_max + 2;
					if (run_parameters.isospin_breaking_1S0==true){
						num_2N_unco_states -= 1;
						num_2N_coup_states  = 1;
					}
				}

				int V_unco_array_size = 0;
				int V_coup_array_size = 0;

				/* Check if we can have coupled channels */
				if (run_parameters.tensor_force){
					V_unco_array_size = Np_WP*Np_WP   * num_2N_unco_states;
					V_coup_array_size = Np_WP*Np_WP*4 * num_2N_coup_states;
				}
				else{
					V_unco_array_size = Np_WP*Np_WP * num_2N_unco_states;
				}

				/* Potential matrices */
				double* V_WP_unco_array = new double [V_unco_array_size];
				double* V_WP_coup_array = new double [V_coup_array_size];

				printf(" - Constructing 2N-potential matrices in WP basis ... \n");
				calculate_potential_matrices_array_in_WP_basis(V_WP_unco_array, num_2N_unco_states,
															V_WP_coup_array, num_2N_coup_states,
															fwp_states,
															pw_states,
															pot_ptr,
															run_parameters);
				printf("   - Done \n");

				/* End of code segment for potential matrix construction */
				/* ################################################################################################################### */
				/* ################################################################################################################### */
				/* ################################################################################################################### */
				/* Start of code segment for scattering wave-packets construction */

				double  E_bound = 0;

				double* e_SWP_unco_array = new double [  (Np_WP+1) * num_2N_unco_states];
				double* e_SWP_coup_array = new double [2*(Np_WP+1) * num_2N_coup_states];

				double* C_WP_unco_array = new double [V_unco_array_size];
				double* C_WP_coup_array = new double [V_coup_array_size];

				swp_statespace swp_states;

				printf(" - Constructing 2N SWPs ... \n");
				make_swp_states(e_SWP_unco_array,
								e_SWP_coup_array,
								C_WP_unco_array,
								C_WP_coup_array,
								V_WP_unco_array,
								V_WP_coup_array,
								num_2N_unco_states,
								num_2N_coup_states,
								E_bound,
								fwp_states,
								swp_states,
								pw_states,
								run_parameters);
				printf("   - Using E_bound = %.5f MeV \n", E_bound);
				printf("   - Done \n");

				/* End of code segment for scattering wave-packets construction */
				/* ################################################################################################################### */
				/* ################################################################################################################### */
				/* ################################################################################################################### */
				/* Start of code segment for locating on-shell nucleon-deuteron states */

				/* Locate and index on-shell bins from input energies */
				find_on_shell_bins(solve_config,
								swp_states,
								run_parameters);

				/* Locate and index deuteron channels */
				find_deuteron_channels(solve_config,
									pw_states);
			
				/* End of code segment for locating on-shell nucleon-deuteron states */
				/* ################################################################################################################### */
				/* ################################################################################################################### */
				/* ################################################################################################################### */

				std::string file_identification =   "_Np_"   + std::to_string(swp_states.Np_WP)
												+ "_Nq_"   + std::to_string(swp_states.Nq_WP)
												+ "_JP_"   + std::to_string(two_J_3N)
												+ "_"      + std::to_string(P_3N)
												+ "_Jmax_" + std::to_string(J_2N_max);
				std::string U_mat_filename_t = run_parameters.output_folder + "/" + "U_PW_elements"
																			+ file_identification
																			+ "_PSI_" + std::to_string(idx_param_set)
																			+ ".txt";

				/* Setup struct containing all indexing relvant for desired on-shell U-matrix solutions */
				channel_os_indexing chn_os_indexing;
				chn_os_indexing.num_T_lab			= solve_config.num_T_lab;
				chn_os_indexing.q_com_idx_array		= solve_config.q_com_idx_array;
				chn_os_indexing.num_deuteron_states = solve_config.deuteron_num_array[chn_3N];
				chn_os_indexing.deuteron_idx_array  = solve_config.deuteron_idx_arrays[chn_3N];

				/* Start of code segment for resolvent matrix (diagonal array) construction */

				/* Resolvent array */
				size_t dense_dim = pw_substates.Nalpha * swp_states.Np_WP * swp_states.Nq_WP;
				cdouble* G_array = new cdouble [dense_dim * chn_os_indexing.num_T_lab];
				
				printf(" - Constructing 3N resolvents ... \n");
				for (int j=0; j<chn_os_indexing.num_T_lab; j++){
					double E = solve_config.E_com_array[j] + swp_states.E_bound;
					calculate_resolvent_array_in_SWP_basis(&G_array[j*dense_dim],
															E,
															swp_states,
															pw_substates,
															run_parameters);
				}
				printf("   - Done \n");
				

				/* End of code segment for resolvent matrix (diagonal array) construction */
				/* ################################################################################################################### */
				/* ################################################################################################################### */
				/* ################################################################################################################### */
				/* Start of code segment for iterations of elastic Faddeev equations */
				
				cdouble* U_array = new cdouble [chn_os_indexing.num_T_lab
											  * chn_os_indexing.num_deuteron_states
											  * chn_os_indexing.num_deuteron_states];

				printf(" - Solving Faddeev equation ... \n");
				solve_faddeev_equations(U_array,
										G_array,
										P123_sparse_val_array,
										P123_sparse_row_array,
										P123_sparse_col_array_csc,
										P123_sparse_dim,
										V_WP_unco_array,
										V_WP_coup_array,
										swp_states,
										chn_os_indexing,
										pw_substates,
										file_identification,
										run_parameters);
				printf("   - Done \n");

				/* End of code segment for iterations of elastic Faddeev equations */
				/* ################################################################################################################### */
				/* ################################################################################################################### */
				/* ################################################################################################################### */
				/* Start of code segment for storing on-shell U-matrix solutions */

				store_U_matrix_elements_txt(U_array,
											solve_config,
											chn_os_indexing,
											run_parameters,
											swp_states,
											pw_substates,
											U_mat_filename_t);

				/* End of code segment for storing on-shell U-matrix solutions */
				/* ################################################################################################################### */
				/* ################################################################################################################### */
				/* ################################################################################################################### */
				/* Delete 3N-channel, model-dependent arrays */
				delete [] V_WP_unco_array;
				delete [] V_WP_coup_array;
				delete [] e_SWP_unco_array;
				delete [] e_SWP_coup_array;
				delete [] C_WP_unco_array;
				delete [] C_WP_coup_array;
				delete [] G_array;
				delete [] U_array;
			}
		}
		/* Delete 3N-channel, model-independent arrays */
		delete [] P123_sparse_val_array;
		delete [] P123_sparse_row_array;
		delete [] P123_sparse_col_array;
	}
	/* End of looping over 3N-channels */

	/* -------------------- End main body of code here -------------------- */

	auto program_end= chrono::system_clock::now();

	chrono::duration<double> total_time = program_end - program_start;
	printf("Total run-time: %.6f\n", total_time.count());

	printf("END OF RUN");
	
	return 0;
}


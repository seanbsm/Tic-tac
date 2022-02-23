
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
#include "make_swp_states.h"
#include "make_resolvent.h"
#include "solve_faddeev.h"
#include "disk_io_routines.h"

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

	/* Current scattering energy */
	size_t  num_T_lab	= 0;
	double* T_lab_array = NULL;
	double* q_com_array = NULL;
	double* E_com_array = NULL;
	
	/* Index lookup arrays for keeping track of on-shell nucleon-deuteron channels in state space */
	int*   q_com_idx_array     = NULL;
	int**  deuteron_idx_arrays = NULL;		// Contains indices of deuteron-channels in given 3N-channel
	int*   deuteron_num_array  = NULL;		// Contains number of deuteron-channels in given 3N-channel

	/* PWE truncation */
	/* Maximum (max) values for J_2N and J_3N (minimum is set to 0 and 1, respectively)*/
	int J_2N_max 	 = run_parameters.J_2N_max;
	int two_J_3N_max = run_parameters.two_J_3N_max;
	if ( two_J_3N_max%2==0 ||  two_J_3N_max<=0 ){
		raise_error("Cannot have even two_J_3N_max");
	}

	/* Wave-packet 3N momentum boundaries */
	double* p_WP_array  = NULL;
	double* q_WP_array  = NULL;

	/* Quadrature 3N momenta per WP cell */
	int Nphi		 = run_parameters.Nphi;
	int Nx 			 = run_parameters.Nx;
	int Np_per_WP	 = run_parameters.Np_per_WP;
	int Nq_per_WP	 = run_parameters.Nq_per_WP;
	double* p_array  = NULL;
	double* q_array  = NULL;
	double* wp_array = NULL;
	double* wq_array = NULL;

	/* Potential matrices */
	double* V_WP_unco_array = NULL;
	double* V_WP_coup_array = NULL;

	/* Quantum numbers of partial-wave expansion in state_3N_array.
	 * All non-specified quantum numbers are either given by nature,
	 * deduced from the values below, or disappear in summations later */
	pw_3N_statespace pw_states;

	/* Three-nucleon channel indexing (a 3N channel is defined by a set {two_J_3N, two_T_3N, P_3N}) */
	int  N_chn_3N 		  = 0;
	int* chn_3N_idx_array = NULL;

	/* Potential model class pointers */
	potential_model* pot_ptr  = potential_model::fetch_potential_ptr(run_parameters.potential_model);

	/* End of code segment for variables and arrays declaration */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* Start of code segment for state space construction */
	
	printf("Constructing 3N partial-wave basis ... \n");
	construct_symmetric_pw_states(N_chn_3N,
								  &chn_3N_idx_array,
								  pw_states,
								  run_parameters);
	printf(" - There are %d 3N-channels \n", N_chn_3N);

	/* Allocate deuteron-channel index-lookup arrays */
	deuteron_idx_arrays = new int* [N_chn_3N];
	deuteron_num_array  = new int  [N_chn_3N];

	/* Small script for finding the largest 3N-channel */
	if (true){
		int largest_Nalpha 	   = 0;
		int largest_Nalpha_idx = 0;
		for (int chn_3N=0; chn_3N<N_chn_3N; chn_3N++){
			int idx_alpha_lower  = chn_3N_idx_array[chn_3N];
			int idx_alpha_upper  = chn_3N_idx_array[chn_3N+1];
			int Nalpha_in_3N_chn = idx_alpha_upper - idx_alpha_lower;

			if (Nalpha_in_3N_chn>largest_Nalpha){
				largest_Nalpha 	   = Nalpha_in_3N_chn;
				largest_Nalpha_idx = chn_3N;
			}
		}
		printf(" - Channel number %d is the largest channel (%d partial-wave states) \n", largest_Nalpha_idx, largest_Nalpha);
	}
	printf(" - Done \n");
	
	printf("Constructing wave-packet (WP) p-momentum bin boundaries ... \n");
	p_WP_array = new double [Np_WP+1];
	make_p_bin_grid(Np_WP, p_WP_array, run_parameters);
	printf(" - Done \n");
	printf("Constructing wave-packet (WP) q-momentum bin boundaries ... \n");
	q_WP_array = new double [Nq_WP+1];
	make_q_bin_grid(Nq_WP, q_WP_array, run_parameters);
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

	printf("Storing q boundaries to CSV-file ... \n");
	//std::string U_mat_foldername = "../../Data/Faddeev_code/U_matrix_elements/";
	std::string q_boundaries_filename = run_parameters.output_folder + "/" + "q_boundaries_Nq_" + to_string(Nq_WP) + ".csv";
	store_q_WP_boundaries_csv(Nq_WP, q_WP_array, q_boundaries_filename);
	printf(" - Done \n");

	/* End of code segment for state space construction */
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

	V_WP_unco_array = new double [V_unco_array_size];
	V_WP_coup_array = new double [V_coup_array_size];
	if (run_parameters.solve_faddeev){
		printf("Constructing 2N-potential matrices in WP basis ... \n");
		calculate_potential_matrices_array_in_WP_basis(V_WP_unco_array, num_2N_unco_states,
													   V_WP_coup_array, num_2N_coup_states,
													   Np_WP, p_WP_array,
													   Np_per_WP, p_array, wp_array,
													   pw_states,
													   pot_ptr,
													   run_parameters);
		printf(" - Done \n");
	}

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
	if (run_parameters.solve_faddeev){
		printf("Constructing 2N SWPs ... \n");
		make_swp_states(e_SWP_unco_array,
						e_SWP_coup_array,
						C_WP_unco_array,
						C_WP_coup_array,
						V_WP_unco_array,
						V_WP_coup_array,
						num_2N_unco_states,
						num_2N_coup_states,
						E_bound,
						Np_WP, p_WP_array,
						pw_states,
						run_parameters);
		printf(" - Using E_bound = %.5f MeV \n", E_bound);
		printf(" - Done \n");
	}

	/* End of code segment for scattering wave-packets construction */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* Start of code segment for storing kinematics of WP statespace */

	double* Eq_WP_boundaries   = new double [Nq_WP+1];
	double* Tlab_WP_boundaries = new double [Nq_WP+1];
	for (size_t q_WP_idx=0; q_WP_idx<Nq_WP+1; q_WP_idx++){
		Eq_WP_boundaries[q_WP_idx]   = com_q_momentum_to_com_energy(q_WP_array[q_WP_idx]);
		Tlab_WP_boundaries[q_WP_idx] = com_momentum_to_lab_energy(q_WP_array[q_WP_idx], E_bound);
	}
	double* q_WP_midpoints     = new double [Nq_WP];
	double* Eq_WP_midpoints    = new double [Nq_WP];
	double* Tlab_WP_midpoints  = new double [Nq_WP];
	for (size_t q_WP_idx=0; q_WP_idx<Nq_WP; q_WP_idx++){
		double Eq_lower = 0.5*(q_WP_array[q_WP_idx]   * q_WP_array[q_WP_idx])  /mu1(E_bound);
		double Eq_upper = 0.5*(q_WP_array[q_WP_idx+1] * q_WP_array[q_WP_idx+1])/mu1(E_bound);
		double E_com = 0.5*(Eq_upper + Eq_lower);
		q_WP_midpoints[q_WP_idx]    = com_energy_to_com_q_momentum(E_com);
		Eq_WP_midpoints[q_WP_idx]   = E_com;
		Tlab_WP_midpoints[q_WP_idx] = com_momentum_to_lab_energy(q_WP_midpoints[q_WP_idx], E_bound);
	}
	printf("Storing kinematic values of WP statespace to txt-file ... \n");
	std::string q_kinematics_filename = run_parameters.output_folder + "/" + "q_kinematics_Nq_" + to_string(Nq_WP) + ".txt";
	store_q_WP_kinematics_txt(Nq_WP,
							  q_WP_array,
							  Eq_WP_boundaries,
							  Tlab_WP_boundaries,
							  q_WP_midpoints,
							  Eq_WP_midpoints,
							  Tlab_WP_midpoints,
							  q_kinematics_filename);
	printf(" - Done \n");

	/* End of code segment for storing kinematics of WP statespace */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* Start of code segment for locating on-shell nucleon-deuteron states */
	if (run_parameters.solve_faddeev){

		/* Use q-momentum bin ENERGY mid-points as on-shell energies if no default input is given */
		std::vector<size_t> q_WP_idx_vec;

		/* Special condition to reduce number of on-shell calculations */
		double Eq_lower = 0;
		double Eq_upper = 0;
		double E_com    = 0;
		double q_m		= 0;
		std::vector<bool>   midpoint_idx_vector   (Nq_WP-1, false);
		std::vector<double> T_lab_midpoint_vector (Nq_WP, false);
		for (size_t q_WP_idx=0; q_WP_idx<Nq_WP; q_WP_idx++){
			Eq_lower = 0.5*(q_WP_array[q_WP_idx]   * q_WP_array[q_WP_idx])  /mu1(E_bound);
			Eq_upper = 0.5*(q_WP_array[q_WP_idx+1] * q_WP_array[q_WP_idx+1])/mu1(E_bound);
			E_com	 = 0.5*(Eq_upper + Eq_lower);
			q_m      = com_energy_to_com_q_momentum(E_com);
			T_lab_midpoint_vector[q_WP_idx] = com_momentum_to_lab_energy(q_m, E_bound);
		}

		double* energy_input_array = NULL;
		int		num_energy_input   = 0;
		read_input_energies(energy_input_array, num_energy_input, run_parameters.energy_input_file);

		//std::vector<double> T_lab_input_list = {1,2,3,4,5,9,10,13,22.7,35,53};
		//std::vector<double> T_lab_input_list = {  3,   4,   5,       6,
		//										  9,  10,  11,      12,
		//										 13,  16,  22.7,    28,
		//										 30,  35,  42,      47.5,
		//										 50,  53,  65,      93.5,
		//										146, 155, 180, 220, 240};
			
		for (size_t q_WP_idx=0; q_WP_idx<Nq_WP-1; q_WP_idx++){
			double T_lab_lower = T_lab_midpoint_vector[q_WP_idx];
			double T_lab_upper = T_lab_midpoint_vector[q_WP_idx+1];
			for (int i=0; i<num_energy_input; i++){
				double T_lab_input = energy_input_array[i];
				/* See if input energy lies between two bin mid-points */
				if (T_lab_lower<=T_lab_input && T_lab_input<=T_lab_upper){
					midpoint_idx_vector[q_WP_idx] = true;
				}
			}
		}
		/* Use on-shell midpoints to set on-shell bins */
		std::vector<bool>   bin_idx_vector   (Nq_WP, false);
		for (size_t q_WP_idx=0; q_WP_idx<Nq_WP-1; q_WP_idx++){
			if (midpoint_idx_vector[q_WP_idx]==true){
				bin_idx_vector[q_WP_idx]   = true;
				bin_idx_vector[q_WP_idx+1] = true;
			}
		}
		/* Append on-shell bin indices to q_WP_idx_vec */
		for (size_t q_WP_idx=0; q_WP_idx<Nq_WP; q_WP_idx++){
			if (bin_idx_vector[q_WP_idx]==true){
				q_WP_idx_vec.push_back(q_WP_idx);
			}
		}
		num_T_lab = q_WP_idx_vec.size();

		T_lab_array = new double [num_T_lab];

		for (size_t Tlab_idx=0; Tlab_idx<num_T_lab; Tlab_idx++){
			size_t q_WP_idx = q_WP_idx_vec[Tlab_idx];
			double Eq_lower = 0.5*(q_WP_array[q_WP_idx]   * q_WP_array[q_WP_idx])  /mu1(E_bound);
			double Eq_upper = 0.5*(q_WP_array[q_WP_idx+1] * q_WP_array[q_WP_idx+1])/mu1(E_bound);
			double E_com 	= 0.5*(Eq_upper + Eq_lower);
			double q 	 	= com_energy_to_com_q_momentum(E_com);
			double T_lab 	= com_momentum_to_lab_energy(q, E_bound);
			T_lab_array[Tlab_idx] = T_lab;
		}

		/* Calculate on-shell momentum and energy in centre-of-mass frame */
		q_com_array = new double [num_T_lab];
		E_com_array = new double [num_T_lab];
		
		for (size_t i=0; i<num_T_lab; i++){
			double q_com = lab_energy_to_com_momentum(T_lab_array[i], E_bound);
			q_com_array[i] = q_com;
			E_com_array[i] = com_q_momentum_to_com_energy(q_com); 
		}

		printf("Locating on-shell q-momentum WP-indices for %zu on-shell energies ... \n", num_T_lab);
		q_com_idx_array = new int [num_T_lab];
		for (size_t idx_Tlab=0; idx_Tlab<num_T_lab; idx_Tlab++){
			double q_com = q_com_array[idx_Tlab];

			int idx_q_bin = -1;
			for (int q_idx_WP=0; q_idx_WP<Nq_WP; q_idx_WP++){
				if (q_WP_array[q_idx_WP]<q_com and q_com<q_WP_array[q_idx_WP+1]){
					idx_q_bin = q_idx_WP;
					break;
				}
			}

			if (idx_q_bin==-1){
				printf("On-shell kinetic energy Tlab=%.3f MeV doesn't exist in WP state space \n", T_lab_array[idx_Tlab]);
				raise_error("Invalid Tlab entered. Exiting ...");
			}
			else{
				q_com_idx_array[idx_Tlab] = idx_q_bin;
			}
		}
		printf(" - On-shell q-momentum WP bins found \n");

		int deuteron_L = 0;
		int deuteron_S = 1;
		int deuteron_J = 1;
		int deuteron_T = 0;
		for (int chn_3N=0; chn_3N<N_chn_3N; chn_3N++){
			/* Lower and upper limits on PW state space for channel */
			int idx_alpha_lower  = chn_3N_idx_array[chn_3N];
			int idx_alpha_upper  = chn_3N_idx_array[chn_3N+1];
			int Nalpha_in_3N_chn = idx_alpha_upper - idx_alpha_lower;

			/* Pointers to sub-arrays of PW state space corresponding to chn_3N */
			int* L_2N_subarray 	   = &pw_states.L_2N_array[idx_alpha_lower];
			int* S_2N_subarray 	   = &pw_states.S_2N_array[idx_alpha_lower];
			int* J_2N_subarray 	   = &pw_states.J_2N_array[idx_alpha_lower];
			int* T_2N_subarray 	   = &pw_states.T_2N_array[idx_alpha_lower];

			std::vector<int> deuteron_chn_indices;

			/* Find indices of deuteron-channels */
			for (int idx_alpha=0; idx_alpha<Nalpha_in_3N_chn; idx_alpha++){
				if (deuteron_L==L_2N_subarray[idx_alpha] &&
				    deuteron_S==S_2N_subarray[idx_alpha] &&
					deuteron_J==J_2N_subarray[idx_alpha] &&
					deuteron_T==T_2N_subarray[idx_alpha]){
					deuteron_chn_indices.push_back(idx_alpha);
				}
			}
			
			/* Copy vector content into array */
			int* chn_3N_deuteron_indices_array = new int [deuteron_chn_indices.size()];
			std::copy(deuteron_chn_indices.begin(), deuteron_chn_indices.end(), chn_3N_deuteron_indices_array);

			/* Add array length and pointer to book-keeping arrays */
			deuteron_idx_arrays[chn_3N] = chn_3N_deuteron_indices_array;
			deuteron_num_array [chn_3N] = deuteron_chn_indices.size();
		}
		printf(" - Nucleon-deuteron channel indices found \n");

		printf(" - Done \n");
	}
	/* End of code segment for locating on-shell nucleon-deuteron states */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* Start of looping over 3N-channels */

	for (int chn_3N=0; chn_3N<N_chn_3N; chn_3N++){

		if (run_parameters.parallel_run==true && chn_3N!=run_parameters.channel_idx){
			continue;
		}

		/* Start of 3N-channel setup */
		/* Lower and upper limits on PW state space for channel */
		int idx_alpha_lower  = chn_3N_idx_array[chn_3N];
		int idx_alpha_upper  = chn_3N_idx_array[chn_3N+1];
		int Nalpha_in_3N_chn = idx_alpha_upper - idx_alpha_lower;

		/* Pointers to sub-arrays of PW state space corresponding to chn_3N */
		pw_3N_statespace pw_substates;
		pw_substates.Nalpha			= idx_alpha_upper - idx_alpha_lower;
		pw_substates.L_2N_array 	= &pw_states.L_2N_array[idx_alpha_lower];
		pw_substates.S_2N_array 	= &pw_states.S_2N_array[idx_alpha_lower];
		pw_substates.J_2N_array 	= &pw_states.J_2N_array[idx_alpha_lower];
		pw_substates.T_2N_array 	= &pw_states.T_2N_array[idx_alpha_lower];
		pw_substates.L_1N_array 	= &pw_states.L_1N_array[idx_alpha_lower];
		pw_substates.two_J_1N_array = &pw_states.two_J_1N_array[idx_alpha_lower];
		pw_substates.two_T_3N_array = &pw_states.two_T_3N_array[idx_alpha_lower];
		pw_substates.two_J_3N_array = &pw_states.two_J_3N_array[idx_alpha_lower];
		pw_substates.P_3N_array		= &pw_states.P_3N_array[idx_alpha_lower];

		/* Channel conserved 3N quantum numbers using first element in channel */
		int two_J_3N = pw_substates.two_J_3N_array[0];
		int P_3N 	 = pw_substates.P_3N_array[0];

		printf("Working on 3N-channel J_3N=%.d/2, PAR=%.d (channel %.d of %.d) with %.d partial-wave states \n", two_J_3N, P_3N, chn_3N+1, N_chn_3N, Nalpha_in_3N_chn);

		/* End of 3N-channel setup */
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
						 Np_WP, p_WP_array,
						 Nq_WP, q_WP_array,
						 Nx,
						 Nphi,
						 J_2N_max,
						 pw_substates,
						 run_parameters,
						 run_parameters.P123_folder);

		/* End of code segment for permutation matrix construction */

		if (run_parameters.solve_faddeev){
			/* Start of code segment for resolvent matrix (diagonal array) construction */

			/* Resolvent array */
			cdouble* G_array = new cdouble [Nalpha_in_3N_chn * Np_WP * Nq_WP * num_T_lab];
			
			printf("Constructing 3N resolvents ... \n");
			for (int j=0; j<num_T_lab; j++){
				double E_com = E_com_array[j];
				//double E_com = 0.03873911528416455;
				double E = E_com + E_bound;
				calculate_resolvent_array_in_SWP_basis(&G_array[j*Nalpha_in_3N_chn*Np_WP*Nq_WP],
													   E,
													   Np_WP,
													   e_SWP_unco_array,
													   e_SWP_coup_array,
													   Nq_WP, q_WP_array,
													   pw_substates,
													   run_parameters);
			}
			printf(" - Done \n");
			

			/* End of code segment for resolvent matrix (diagonal array) construction */
			/* Start of code segment for iterations of elastic Faddeev equations */
			int 	 num_deuteron_states = deuteron_num_array[chn_3N];
			int* 	 deuteron_idx_array  = deuteron_idx_arrays[chn_3N];
			
			cdouble* U_array = new cdouble [num_T_lab * num_deuteron_states * num_deuteron_states];

            std::string file_identification =   "_Np_" + std::to_string(Np_WP)
									          + "_Nq_"   + std::to_string(Nq_WP)
									          + "_JP_"   + std::to_string(two_J_3N)
									          + "_"      + std::to_string(P_3N)
									          + "_Jmax_" + std::to_string(J_2N_max);

			printf("Solving Faddeev equations ... \n");
			solve_faddeev_equations(U_array,
									G_array,
									P123_sparse_val_array,
									P123_sparse_row_array,
									P123_sparse_col_array,
									P123_sparse_dim,
									C_WP_unco_array,
									C_WP_coup_array,
									V_WP_unco_array,
									V_WP_coup_array,
									num_2N_unco_states,
									num_2N_coup_states,
									q_com_idx_array, num_T_lab,
									deuteron_idx_array, num_deuteron_states,
									Nq_WP,
									Np_WP,
									pw_substates,
									file_identification,
					                run_parameters);
			printf(" - Done \n");

			/* End of code segment for iterations of elastic Faddeev equations */
			/* Start of code segment for storing on-shell U-matrix solutions */

			//std::string U_mat_filename = run_parameters.output_folder + "/" + "U_PW_elements"
			//                                                                + file_identification
			//														        + ".csv";
			//store_U_matrix_elements_csv(U_array,
			//						    q_com_idx_array,    (size_t) num_T_lab,
			//				  		    deuteron_idx_array, (size_t) num_deuteron_states,
			//						    pw_substates,
			//						    U_mat_filename);
			
			std::string U_mat_filename_t = run_parameters.output_folder + "/" + "U_PW_elements"
			                                                                  + file_identification
																	          + ".txt";
			store_U_matrix_elements_txt(U_array,
										run_parameters.potential_model,
										Np_WP,
										Nq_WP,
										E_bound,
										T_lab_array,
										E_com_array,
									    q_com_idx_array,    (size_t) num_T_lab,
							  		    deuteron_idx_array, (size_t) num_deuteron_states,
									    pw_substates,
									    U_mat_filename_t);

			/* End of code segment for storing on-shell U-matrix solutions */
		}
	}
	/* End of looping over 3N-channels */

	/* -------------------- End main body of code here -------------------- */

	auto program_end= chrono::system_clock::now();

	chrono::duration<double> total_time = program_end - program_start;
	printf("Total run-time: %.6f\n", total_time.count());

	printf("END OF RUN");
	
	return 0;
}


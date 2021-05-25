
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
#include "General_functions/kinetic_conversion.h"
#include "Interactions/potential_model.h"
#include "make_potential_matrix.h"
#include "make_swp_states.h"
#include "make_resolvent.h"
#include "solve_faddeev.h"
#include "basis_transformations.h"
#include "store_functionality.h"

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
 
 - Call function solve_faddeev_equation (Nalpha_chn = number of alpha in current 3N channel)
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
	
	/* Array-job information for simple parallellism on several nodes */
	int  job_ID       = 0;
	int  num_jobs     = 0;
	bool job_array_on = false;
	if (argc==3){
		job_ID       = atoi(argv[1]);
		num_jobs     = atoi(argv[2]);
		job_array_on = true;
		if (job_ID >= num_jobs){
			raise_error("Job_ID is greater than or equal to the number of jobs");
		}
	}
	else if (argc!=1){
		raise_error("Invalid number of input arguments in main.cpp");
	}

	auto program_start = chrono::system_clock::now();
	
	/* ------------------- Start main body of code here ------------------- */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* Start of code segment for parameters, variables and arrays declaration */

	bool default_Tlab_input = false;
	double Tlab_max = 100;

	/* Current scattering energy */
	size_t  num_T_lab	= 0;
	double* T_lab_array = NULL;//new double [5] {7.49793, 12.7426, 22.4742, 37.4248, 52.3506};
	double* q_com_array = NULL;
	double* E_com_array = NULL;
	
	/* Use Tlab = 1, 2, 3, ..., num_T_lab (TEMPORARY IMPLEMENTATION) */
	if (default_Tlab_input){
		num_T_lab = 1;
		T_lab_array = new double [num_T_lab];
		for (size_t i=0; i<num_T_lab; i++){
			T_lab_array[i] = i + 0.0001*(i==0);
			//T_lab_array[i] = 3*i + 0.0001*(i==0);
		}
	}
	
	/* Index lookup arrays for keeping track of on-shell nucleon-deuteron channels in state space */
	int*   q_com_idx_array     = NULL;
	int**  deuteron_idx_arrays = NULL;		// Contains indices of deuteron-channels in given 3N-channel
	int*   deuteron_num_array  = NULL;		// Contains number of deuteron-channels in given 3N-channel

	/* Setting to store calculated P123 matrix in WP basis to h5-file */
	bool calculate_and_store_P123 = false;
	/* Setting to solve Faddeev or not. Handy if we only want to
	 * precalculate permutation matrices, or to calculate both permutation matrices
	 * and solve Faddeev in a single run */
	bool solve_faddeev		      = true;
	bool production_run			  = true;

	/* PWE truncation */
	/* Maximum (max) values for J_2N and J_3N (minimum is set to 0 and 1, respectively)*/
	int J_2N_max 	 = 1;
	int two_J_3N_max = 3;
	if ( two_J_3N_max%2==0 ||  two_J_3N_max<=0 ){
		raise_error("Cannot have even two_J_3N_max");
	}

	/* Wave-packet 3N momenta */
	int Np_WP	   	 = 8;
	int Nq_WP	   	 = 8;
	double* p_WP_array  = NULL;
	double* q_WP_array  = NULL;

	/* Quadrature 3N momenta per WP cell */
	int Nphi		 = 50;
	int Nx 			 = 20;
	int Np_per_WP	 = 8;
	int Nq_per_WP	 = 8;
	double* p_array  = NULL;
	double* q_array  = NULL;
	double* wp_array = NULL;
	double* wq_array = NULL;

	/* Potential matrices */
	double* V_WP_unco_array = NULL;
	double* V_WP_coup_array = NULL;
	double** V_WP_2N_arrays = NULL;

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

	/* End of code segment for variables and arrays declaration */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
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
	make_p_bin_grid(Np_WP, p_WP_array);
	printf(" - Done \n");
	printf("Constructing wave-packet (WP) q-momentum bin boundaries ... \n");
	q_WP_array = new double [Nq_WP+1];
	make_q_bin_grid(Nq_WP, q_WP_array);
	printf(" - Done \n");

	//for (int i=0; i<Nq_WP+1; i++){
	//	printf("%.3f\n", q_WP_array[i]);
	//}
	//return 0;

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
	std::string U_mat_foldername = "../../Data/U_matrix_elements/";
	std::string q_boundaries_filename = U_mat_foldername + "q_boundaries_Nq_" + to_string(Nq_WP) + ".csv";
	store_q_WP_boundaries_csv(Nq_WP, q_WP_array,
						      q_boundaries_filename);
	printf(" - Done \n");

	/* End of code segment for state space construction */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* Start of code segment for locating on-shell nucleon-deuteron states */
	if (solve_faddeev){

		/* Use q-momentum bin mid-points as on-shell energies if no default input is given */
		if (default_Tlab_input==false){
			for (size_t q_WP_idx=0; q_WP_idx<Nq_WP; q_WP_idx++){
				double q_upper_boundary = q_WP_array[q_WP_idx+1];
				double T_lab = com_momentum_to_lab_energy(q_upper_boundary);
				if (T_lab>Tlab_max){
					break;
				}
				else{
					num_T_lab += 1;
				}
			}

			T_lab_array = new double [num_T_lab];

			for (size_t q_WP_idx=0; q_WP_idx<num_T_lab; q_WP_idx++){
				double q_mid_point = 0.5*(q_WP_array[q_WP_idx+1] + q_WP_array[q_WP_idx]);
				double T_lab = com_momentum_to_lab_energy(q_mid_point);
				
				T_lab_array[q_WP_idx] = T_lab;
			}
		}

		/* Calculate on-shell momentum and energy in centre-of-mass frame */
		//double mu	 = 2*Mp*Md/(Mp+Md);
		q_com_array = new double [num_T_lab];
		E_com_array = new double [num_T_lab];
		//double mu	 = 2*Mp*Md/(Mp+Md);
		for (size_t i=0; i<num_T_lab; i++){
			double q_com = lab_energy_to_com_momentum(T_lab_array[i]);
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
			int* L_2N_subarray 	   = &L_2N_array[idx_alpha_lower];
			int* S_2N_subarray 	   = &S_2N_array[idx_alpha_lower];
			int* J_2N_subarray 	   = &J_2N_array[idx_alpha_lower];
			int* T_2N_subarray 	   = &T_2N_array[idx_alpha_lower];

			std::vector<int> deuteron_chn_indices;

			/* Find indices of deuteron-channels */
			for (int idx_alpha=0; idx_alpha<Nalpha_in_3N_chn; idx_alpha++){
				if (deuteron_L==L_2N_subarray[idx_alpha] &
				    deuteron_S==S_2N_subarray[idx_alpha] &
					deuteron_J==J_2N_subarray[idx_alpha] &
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
	/* Start of code segment for potential matrix construction */

	int V_unco_array_size = Np_WP*Np_WP   * 2*(J_2N_max+1);
	int V_coup_array_size = Np_WP*Np_WP*4 *    J_2N_max;
	V_WP_unco_array = new double [V_unco_array_size];
	V_WP_coup_array = new double [V_coup_array_size];
	V_WP_2N_arrays  = new double* [Nalpha * Nalpha];
	if (solve_faddeev){
		//pot_ptr_np = potential_model::fetch_potential_ptr("LO_internal", "np");
		//pot_ptr_nn = potential_model::fetch_potential_ptr("LO_internal", "nn");
		//pot_ptr_np = potential_model::fetch_potential_ptr("N2LOopt", "np");
		//pot_ptr_nn = potential_model::fetch_potential_ptr("N2LOopt", "nn");
		//pot_ptr_np = potential_model::fetch_potential_ptr("Idaho_N3LO", "np");
		//pot_ptr_nn = potential_model::fetch_potential_ptr("Idaho_N3LO", "nn");
		pot_ptr_np = potential_model::fetch_potential_ptr("malfliet_tjon", "np");
		pot_ptr_nn = potential_model::fetch_potential_ptr("malfliet_tjon", "nn");
	
		//double temparray [6];
		//double qi = 5; double qo=10;
		//int L=0; int S3=1; int J3=1; int S1=0; int J1=0;
		//std::cout << "Potential vals:" << std::endl;
		//pot_ptr_nn->V(qi, qo, false, L,S1,J1, temparray);
		//std::cout << temparray[0] << std::endl;
		//pot_ptr_nn->V(qi, qo, true,  L,S3,J3, temparray);
		//std::cout << temparray[2] << std::endl;
		//return 0;

		//for (int i=0; i<Np_WP+1; i++){
		//	std::cout << p_WP_array[i] << std::endl;
		//}
		//return 0;

		printf("Constructing 2N-potential matrices in WP basis ... \n");
		calculate_potential_matrices_array_in_WP_basis(V_WP_unco_array,
													   V_WP_coup_array,
													   V_WP_2N_arrays,
													   true,
													   Np_WP, p_WP_array,
													   Np_per_WP, p_array, wp_array,
													   Nalpha, L_2N_array, S_2N_array, J_2N_array, T_2N_array,
													   J_2N_max,
													   pot_ptr_nn,
													   pot_ptr_np);
		printf(" - Done \n");
	}

	/* End of code segment for potential matrix construction */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* Start of code segment for scattering wave-packets construction */

	double  E_bound = 0;

	double* e_SWP_unco_array = new double [  (Np_WP+1) * 2*(J_2N_max+1)];
	double* e_SWP_coup_array = new double [2*(Np_WP+1) *    J_2N_max];

	double* C_WP_unco_array = new double [V_unco_array_size];
	double* C_WP_coup_array = new double [V_coup_array_size];
	if (solve_faddeev){
		printf("Constructing 2N SWPs ... \n");
		make_swp_states(e_SWP_unco_array,
						e_SWP_coup_array,
						C_WP_unco_array,
						C_WP_coup_array,
						V_WP_unco_array,
						V_WP_coup_array,
						E_bound,
						Np_WP, p_WP_array,
						Nalpha, L_2N_array, S_2N_array, J_2N_array, T_2N_array,
						J_2N_max);
		printf(" - Using E_bound = %.5f MeV \n", E_bound);
		printf(" - Done \n");

		for (int i=0; i<Np_WP+1; i++){
			printf("%.3f\n", p_WP_array[i]);
		}
		printf("\n");
		for (int i=0; i<2*(Nq_WP+1); i++){
			printf("%.3f\n", e_SWP_coup_array[i]);
		}
		return 0;
	}

	//for (int i=2; i<Nalpha; i++){
	//	for (int j=0; j<2*Np_WP; j++){
	//		printf("%.3f\n", e_SWP_unco_array[j]);
	//	}
	//}

	/* End of code segment for scattering wave-packets construction */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* Start of looping over 3N-channels */

	for (int chn_3N=0; chn_3N<N_chn_3N; chn_3N++){

		/* Check channel distribution in case of parallell execution */
		if (job_array_on){
			int num_chn_per_node = N_chn_3N/num_jobs;

			int chn_lower = job_ID*num_chn_per_node;
			int chn_upper = (job_ID+1)*num_chn_per_node;

			/* Last job handles the last number of permutation matrices - NOT OPTIMAL*/
			if (job_ID == num_jobs-1){
				chn_upper = N_chn_3N;
			}
			
			/* Skip to next chn-iteration if chn does not fit in current job_ID's domain */
			if (chn_3N<chn_lower || chn_upper<=chn_3N){
				continue;
			}
		}

		/* Start of 3N-channel setup */
		/* Lower and upper limits on PW state space for channel */
		int idx_alpha_lower  = chn_3N_idx_array[chn_3N];
		int idx_alpha_upper  = chn_3N_idx_array[chn_3N+1];
		int Nalpha_in_3N_chn = idx_alpha_upper - idx_alpha_lower;

		/* Channel conserved 3N quantum numbers using first element in channel */
		int two_J_3N = two_J_3N_array[idx_alpha_lower];
		int two_T_3N = two_T_3N_array[idx_alpha_lower];
		int P_3N 	 = P_3N_array    [idx_alpha_lower];

		/* Pointers to sub-arrays of PW state space corresponding to chn_3N */
		int* L_2N_subarray 	   = &L_2N_array[idx_alpha_lower];
		int* S_2N_subarray 	   = &S_2N_array[idx_alpha_lower];
		int* J_2N_subarray 	   = &J_2N_array[idx_alpha_lower];
		int* T_2N_subarray 	   = &T_2N_array[idx_alpha_lower];
		int* L_1N_subarray 	   = &L_1N_array[idx_alpha_lower];
		int* two_J_1N_subarray = &two_J_1N_array[idx_alpha_lower];

		printf("Working on 3N-channel J_3N=%.d/2, T_3N=%.d/2, PAR=%.d (channel %.d of %.d) with %.d partial-wave states \n", two_J_3N, two_T_3N, P_3N, chn_3N+1, N_chn_3N, Nalpha_in_3N_chn);

		/* End of 3N-channel setup */
		/* Start of code segment for permutation matrix construction */
		double* P123_sparse_val_array = NULL;
		int* 	P123_sparse_row_array = NULL;
		int* 	P123_sparse_col_array = NULL;
		size_t	P123_sparse_dim		  = 0;

		/* Default filename for current chn_3N - used for storage and reading P123 */
		
		std::string permutation_matrices_folder = "../../Data/permutation_matrices/";
		std::string P123_filename =    permutation_matrices_folder + "P123_sparse_JTP_"
									 + to_string(two_J_3N) + "_" + to_string(two_T_3N) + "_" + to_string(P_3N)
									 + "_Np_" + to_string(Np_WP) + "_Nq_" + to_string(Nq_WP)
									 + "_J2max_" + to_string(J_2N_max) + ".h5";
		if (chn_3N!=2){
			continue;
		}
		if (calculate_and_store_P123){
			double* x_array  = new double [Nx];
			double* wx_array = new double [Nx];
			gauss(x_array, wx_array, Nx);
	
			printf("Calculating P123 ... \n");
			auto timestamp_P123_calc_start = chrono::system_clock::now();
			calculate_permutation_matrices_for_all_3N_channels(&P123_sparse_val_array,
															   &P123_sparse_row_array,
															   &P123_sparse_col_array,
															   P123_sparse_dim,
															   production_run,
															   Np_WP, p_WP_array,
														   	   Nq_WP, q_WP_array,
														   	   Nx, x_array, wx_array,
															   Nphi,
															   J_2N_max,
														   	   Nalpha_in_3N_chn,
														   	   L_2N_subarray,
							  								   S_2N_subarray,
							  								   J_2N_subarray,
							  								   T_2N_subarray,
							  								   L_1N_subarray,
							  								   two_J_1N_subarray,
														 	   two_J_3N,
														 	   two_T_3N,
															   P_3N);
			auto timestamp_P123_calc_end = chrono::system_clock::now();
			chrono::duration<double> time_P123_calc = timestamp_P123_calc_end - timestamp_P123_calc_start;
			printf(" - Done. Time used: %.6f\n", time_P123_calc.count());
	
			printf("Storing P123 to h5 ... \n");
			auto timestamp_P123_store_start = chrono::system_clock::now();
			store_sparse_permutation_matrix_for_3N_channel_h5(P123_sparse_val_array,
															  P123_sparse_row_array,
															  P123_sparse_col_array,
															  P123_sparse_dim,
															  Np_WP, p_WP_array,
															  Nq_WP, q_WP_array,
															  Nalpha_in_3N_chn,
															  L_2N_subarray,
															  S_2N_subarray,
															  J_2N_subarray,
															  T_2N_subarray,
															  L_1N_subarray,
															  two_J_1N_subarray,
															  two_J_3N,
															  two_T_3N,
															  P_3N,
													   		  P123_filename,
															  true);
			auto timestamp_P123_store_end = chrono::system_clock::now();
			chrono::duration<double> time_P123_store = timestamp_P123_store_end - timestamp_P123_store_start;
			printf(" - Done. Time used: %.6f\n", time_P123_store.count());
		}
		else if (solve_faddeev){
			//std::string P123_filename_2 =    "../../Data/permutation_matrices/P123_sparse_JTP_"
			//						 + to_string(two_J_3N) + "_" + to_string(two_T_3N) + "_" + to_string(P_3N)
			//						 + "_Np_" + to_string(Np_WP) + "_Nq_" + to_string(Nq_WP)
			//						 + "_J2max_" + to_string(J_2N_max) + ".h5";
			printf("Reading P123 from h5 ... \n");
			
			//double* P123_sparse_val_array_t = NULL;
			//int* 	P123_sparse_row_array_t = NULL;
			//int* 	P123_sparse_col_array_t = NULL;
			//size_t  P123_sparse_dim_t		= 0;

			auto timestamp_P123_read_start = chrono::system_clock::now();
			read_sparse_permutation_matrix_for_3N_channel_h5( &P123_sparse_val_array,
															  &P123_sparse_row_array,
															  &P123_sparse_col_array,
															  P123_sparse_dim,
															  Np_WP, p_WP_array,
															  Nq_WP, q_WP_array,
															  Nalpha_in_3N_chn,
															  L_2N_subarray,
															  S_2N_subarray,
															  J_2N_subarray,
															  T_2N_subarray,
															  L_1N_subarray,
															  two_J_1N_subarray,
															  two_J_3N,
															  two_T_3N,
															  P_3N,
													   		  P123_filename,
															  true);
			auto timestamp_P123_read_end = chrono::system_clock::now();
			chrono::duration<double> time_P123_read = timestamp_P123_read_end - timestamp_P123_read_start;
			printf(" - Done. Time used: %.6f\n", time_P123_read.count());
			
			//if (P123_sparse_dim_t==P123_sparse_dim){
			//	//int row_idx = 0;
			//	for (int idx=0; idx<P123_sparse_dim; idx++){
			//		//if (P123_sparse_row_array[idx]==row_idx){
			//		bool check1 = (abs(P123_sparse_val_array_t[idx]-P123_sparse_val_array[idx])>1e-15);
			//		bool check2 = (P123_sparse_row_array_t[idx]!=P123_sparse_row_array[idx]);
			//		bool check3 = (P123_sparse_col_array_t[idx]!=P123_sparse_col_array[idx]);
			//		if (check1||check2||check3){
			//			std::cout << "Value wrong, idx: " << idx << std::endl;
			//			std::cout << "BM val:   " << P123_sparse_val_array_t[idx] << std::endl;
			//			std::cout << "BM row:   " << P123_sparse_row_array_t[idx] << std::endl;
			//			std::cout << "BM col:   " << P123_sparse_col_array_t[idx] << std::endl;
			//			std::cout << "Prog val: " << P123_sparse_val_array[idx] << std::endl;
			//			std::cout << "Prog row: " << P123_sparse_row_array[idx] << std::endl;
			//			std::cout << "Prog col: " << P123_sparse_col_array[idx] << std::endl;
			//			raise_error("element mismatch");
			//		}
			//		//}
			//	}
			//}
			//else{
			//	std::cout << "BM dim:   " << P123_sparse_dim_t << std::endl;
			//	std::cout << "Prog dim: " << P123_sparse_dim << std::endl;
			//	raise_error("dim not right");
			//}
		}
		/* End of code segment for permutation matrix construction */

		if (solve_faddeev){
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
													   Nalpha_in_3N_chn,
													   L_2N_subarray,
													   S_2N_subarray,
													   J_2N_subarray,
													   T_2N_subarray);
				//printf("%.10f\n\n",com_energy_to_com_q_momentum(3.8306227657971945-1.9624153103111233, -1.9624153103111233));
				//double q_com = com_energy_to_com_q_momentum(E_com);
				//printf("Ecm = %.15f MeV | q = %.15f MeV:\n", E_com, q_com);
				//for (int i=0; i<Np_WP*Nq_WP; i++){
				//	printf("%.15e %.15e\n", G_array[i]);
				//}
			}
			//return 0;
			printf(" - Done \n");
			

			/* End of code segment for resolvent matrix (diagonal array) construction */
			/* Start of code segment for iterations of elastic Faddeev equations */
			int 	 num_deuteron_states = deuteron_num_array[chn_3N];
			int* 	 deuteron_idx_array  = deuteron_idx_arrays[chn_3N];
			
			cdouble* U_array = new cdouble [num_T_lab * num_deuteron_states * num_deuteron_states];

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
									q_com_idx_array, num_T_lab,
									deuteron_idx_array, num_deuteron_states,
									J_2N_max,
									Nq_WP,
									Np_WP,
									Nalpha_in_3N_chn,
									L_2N_subarray,
									S_2N_subarray,
									J_2N_subarray,
									T_2N_subarray,
									L_1N_subarray, 
									two_J_1N_subarray);
			printf(" - Done \n");

			/* End of code segment for iterations of elastic Faddeev equations */
			/* Start of code segment for storing on-shell U-matrix solutions */



			/* End of code segment for storing on-shell U-matrix solutions */
		}
	}
	/* End of looping over 3N-channels */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* ################################################################################################################### */
	/* Start of code segment for calculating scattering observables */
	
	/* End of code segment for calculating scattering observables */


	/* -------------------- End main body of code here -------------------- */

	auto program_end= chrono::system_clock::now();

	chrono::duration<double> total_time = program_end - program_start;
	printf("Total run-time: %.6f\n", total_time.count());

	printf("END OF RUN");
	
	return 0;
}



#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

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

using namespace std;

int main(int argc, char* argv[]){

	auto program_start = chrono::system_clock::now();
	
	/* ------------------- Start main body of code here ------------------- */
	/* Start of code segment for parameters, variables and arrays declaration */

	/* Current scattering energy */
	double E = 10;

	/* Tritium bound-state quantum numbers */
	int two_J_3N  	 = 1;
    int two_T_3N  	 = 1;
    int parity_3N 	 = 1;

	/* PWE truncation */
	/* Maximum (max) and minimum (min) values for J_2N and J_1N */
    int J_2N_min 	 = 0;	// The LS-solver will fail if this is not zero - I haven't taken this into account in my indexing
    int J_2N_max 	 = 3;

	/* Wave-packet 3N momenta */
	int Np_WP	   	 = 96;
	int Nq_WP	   	 = 6;
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

	/* Quantum numbers of partial-wave expansion in state_3N_array */
	int  Nalpha   		= 0;		  // Number of partial waves, given by construct_symmetric_pw_states
    int* L_2N_array     = NULL;     // pair angular momentum
    int* S_2N_array     = NULL;     // pair total spin
    int* J_2N_array     = NULL;     // pair total angular momentum
    int* T_2N_array     = NULL;     // pair total isospin
    int* l_3N_array     = NULL;     // three-nucleon angular momentum (?)
    int* two_j_3N_array = NULL; 	  // three-nucleon total angular momentum x2 (?)

	/* Potential model class pointers */
	potential_model* pot_ptr_np = NULL;
	potential_model* pot_ptr_nn = NULL;

	/* End of code segment for variables and arrays decleration */
	/* Start of code segment for state space construction */

	cout << "Constructing 3N partial-wave basis" << endl;
	construct_symmetric_pw_states(two_J_3N, two_T_3N, parity_3N,
								  J_2N_min, J_2N_max,
								  Nalpha, &L_2N_array, &S_2N_array, &J_2N_array, &T_2N_array, &l_3N_array, &two_j_3N_array);

	cout << "Constructing wave-packet (WP) p-momentum bin boundaries" << endl;
	p_WP_array = new double [Np_WP+1];
	make_p_bin_grid(Np_WP, p_WP_array);
	cout << "Constructing wave-packet (WP) q-momentum bin boundaries" << endl;
	q_WP_array = new double [Nq_WP+1];
	make_q_bin_grid(Nq_WP, q_WP_array);

	cout << "Constructing p quadrature mesh per WP, for all WPs" << endl;
	p_array  = new double [Np_per_WP*Np_WP];
	wp_array = new double [Np_per_WP*Np_WP];
	make_p_bin_quadrature_grids(Np_WP, p_WP_array,
                                Np_per_WP, p_array, wp_array);
	cout << "Constructing q quadrature mesh per WP, for all WPs" << endl;
	q_array  = new double [Nq_per_WP*Nq_WP];
	wq_array = new double [Nq_per_WP*Nq_WP];
	make_q_bin_quadrature_grids(Nq_WP, q_WP_array,
                                Nq_per_WP, q_array, wq_array);

	/* End of code segment for state space construction */
	/* Start of code segment for permutation matrix construction */

	cout << "Read P123 dimensions from h5-file" << endl;
	bool interpolate_P123_matrix = false;
	if (interpolate_P123_matrix){
		int Nalpha_P123 = 0;
		int Np_P123 = 0;
		int Nq_P123 = 0;
		get_h5_P123_dimensions(P123_file_path, Nalpha_P123, Np_P123, Nq_P123);
    	P123_array 	  = new double [Np_P123 * Nq_P123 * Nalpha_P123 * Np_P123 * Nq_P123 * Nalpha_P123];
		P123_p_array  = new double [Np_P123];
		P123_q_array  = new double [Nq_P123];
		P123_wp_array = new double [Np_P123];
		P123_wq_array = new double [Nq_P123];
		cout << "Reading P123 from file" << endl;
		read_P123_h5_data_file(P123_file_path,
							   P123_array,
							   Np_P123, P123_p_array, P123_wp_array,
							   Nq_P123, P123_q_array, P123_wq_array);
	}
	else{
		int Nx 			 = 20;
		double* x_array  = new double [Nx];
    	double* wx_array = new double [Nx];
		gauss(x_array, wx_array, Nx);

		P123_array = new double [Nalpha*Np_WP*Nq_WP * Nalpha*Np_WP*Nq_WP];

		//calculate_permutation_matrix(P123_array,
        //                          	 Nq_WP*Nq_per_WP, q_array, wq_array, Np_per_WP, Np_WP, p_WP_array,
        //                          	 Np_WP*Np_per_WP, p_array, wp_array, Nq_per_WP, Nq_WP, q_WP_array,
        //                          	 Nx, x_array, wx_array,
        //                          	 Nalpha,
        //                          	 L_2N_array,
        //                          	 S_2N_array,
        //                          	 J_2N_array,
        //                          	 T_2N_array,
        //                          	 l_3N_array,
        //                          	 two_j_3N_array,
        //                          	 two_J_3N, two_T_3N, parity_3N);
	}

	/* End of code segment for permutation matrix construction */
	/* Start of code segment for potential matrix construction */

	//pot_ptr_np = potential_model::fetch_potential_ptr("LO_internal", "np");
	//pot_ptr_nn = potential_model::fetch_potential_ptr("LO_internal", "nn");
	//pot_ptr_np = potential_model::fetch_potential_ptr("N2LOopt", "np");
	//pot_ptr_nn = potential_model::fetch_potential_ptr("N2LOopt", "nn");
	pot_ptr_np = potential_model::fetch_potential_ptr("Idaho_N3LO", "np");
	pot_ptr_nn = potential_model::fetch_potential_ptr("Idaho_N3LO", "nn");

	cout << "Constructing 2N-potential matrices in WP basis" << endl;
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

	/* End of code segment for potential matrix construction */
	/* Start of code segment for scattering wave-packets construction */

	double* p_SWP_unco_array = new double [  (Np_WP+1) * 2*(J_2N_max+1)];
	double* p_SWP_coup_array = new double [2*(Np_WP+1) *    J_2N_max];

	double* C_WP_unco_array = new double [V_unco_array_size];
	double* C_WP_coup_array = new double [V_coup_array_size];
	
	cout << "Constructing 2N SWPs" << endl;
	make_swp_states(p_SWP_unco_array,
					p_SWP_coup_array,
					C_WP_unco_array,
					C_WP_coup_array,
					V_WP_unco_array,
                    V_WP_coup_array,
					Np_WP, p_WP_array,
					Nalpha, L_2N_array, S_2N_array, J_2N_array, T_2N_array,
					J_2N_max);

	/* End of code segment for scattering wave-packets construction */
	/* Start of code segment for resolvent matrix (diagonal array) construction */

	/* Resolvent array */
	cdouble* G_array = new cdouble [Nalpha * Nq_WP * Np_WP];

	cout << "Constructing 3N resolvents" << endl;
	calculate_resolvent_array_in_SWP_basis(G_array,
                                           E,
                                           Np_WP,
                                           p_SWP_unco_array,
					                       p_SWP_coup_array,
					                       Nq_WP, q_WP_array,
					                       Nalpha, L_2N_array, S_2N_array, J_2N_array, T_2N_array);

	/* End of code segment for resolvent matrix (diagonal array) construction */
	/* Start of code segment for iterations of elastic Faddeev equations */

	cdouble* U_array = new cdouble [Nalpha*Np_WP*Nq_WP * Nalpha*Np_WP*Nq_WP];

	/* End of code segment for iterations of elastic Faddeev equations */
	/* Start of code segment for calculating scattering observables */


	/* End of code segment for calculating scattering observables */


	/* -------------------- End main body of code here -------------------- */

	auto program_end= chrono::system_clock::now();

	chrono::duration<double> total_time = program_end - program_start;
	cout << "Total run-time:                         " << total_time.count() << endl;

	std::cout << "END OF RUN" << std::endl;
	
    return 0;
}


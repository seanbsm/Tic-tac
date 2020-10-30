
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

/* Time-keeping modules */
#include <chrono>
#include <ctime>

#include "make_pw_symm_states.h"
#include "kinetics.h"
#include "permutation_operators.h"
#include "state_antisymmetrization.h"
#include "faddeev_iterator.h"
#include "General_functions/gauss_legendre.h"
#include "Interactions/potential_model.h"
#include "Triton_states/read_psi.h"

#include "lippmann_schwinger_solver.h"

using namespace std;

int main(int argc, char* argv[]){

	auto program_start = chrono::system_clock::now();
	
	/* Start main body of code here */

	/* Tritium bound-state quantum numbers */
	int two_J_3N  	 = 1;
    int two_T_3N  	 = 1;
    int parity_3N 	 = 1;

	/* PWE truncation */
	/* Maximum (max) and minimum (min) values for J_2N and J_1N */
    int two_J_1N_min = 1;
    int two_J_1N_max = 4;
    int J_2N_min 	 = 0;	// The LS-solver will fail if this is not zero - I haven't taken this into account in my indexing
    int J_2N_max 	 = 3;

	/* Quadrature 3N momenta */
	int Np	   		 = 32;
	int Nq	   		 = 10;
	int Nx 			 = 20;
	int Nalpha 		 =  0;
	double* p_array  = NULL;
	double* q_array  = NULL;
	double* wp_array = NULL;
	double* wq_array = NULL;
    double* x_array  = new double [Nx];
    double* wx_array = new double [Nx];
	calculate_angular_quadrature_grids(x_array, wx_array, Nx);

	/* Momentum-representation of 3N state at quadrature nodes in p_array and q_array */
	double* state_3N_symm_array = NULL;
	double* state_3N_asym_array = NULL;

	/* Permutation operator */
	double* P123_array  = NULL;

	/* Anti-symmetrization operator */
	double* A123  = NULL;

	/* Quantum numbers of partial-wave expansion in state_3N_array */
    int* L_2N     = NULL;     // pair angular momentum
    int* S_2N     = NULL;     // pair total spin
    int* J_2N     = NULL;     // pair total angular momentum
    int* T_2N     = NULL;     // pair total isospin
    int* l_3N     = NULL;     // three-nucleon angular momentum (?)
    int* two_j_3N = NULL; 	  // three-nucleon total angular momentum x2 (?)

	/* Tells the program to read pre-calculated antisymmetric triton states.
	 * Handy for small tests since the P123-file can be huge */
	bool use_premade_symmetric_states 	  = false;
	bool use_premade_antisymmetric_states = false;

	//potential_model* pot_ptr_np = potential_model::fetch_potential_ptr("LO_internal", "np");
	//potential_model* pot_ptr_nn = potential_model::fetch_potential_ptr("LO_internal", "nn");
	potential_model* pot_ptr_np = potential_model::fetch_potential_ptr("Idaho_N3LO", "np");
	potential_model* pot_ptr_nn = potential_model::fetch_potential_ptr("Idaho_N3LO", "nn");
	
	if (use_premade_symmetric_states){
		cout << "Reading states from file" << endl;
		get_all_states(&state_3N_symm_array, &state_3N_asym_array, Np, &p_array, &wp_array, Nq, &q_array, &wq_array, Nalpha, &L_2N, &S_2N, &J_2N, &T_2N, &l_3N, &two_j_3N);
	}
	else{
		cout << "Constructing 3N partial-wave basis" << endl;
		construct_symmetric_pw_states(two_J_3N, two_T_3N, parity_3N,
									  two_J_1N_min, two_J_1N_max, J_2N_min, J_2N_max,
									  Nalpha, &L_2N, &S_2N, &J_2N, &T_2N, &l_3N, &two_j_3N);
		cout << "Constructing p mesh" << endl;
		p_array  = new double [Np];
		wp_array = new double [Np];
		gauss(p_array, wp_array, Np); rangeChange_0_inf(p_array, wp_array, 1000., Np);
		cout << "Constructing q mesh" << endl;
		q_array  = new double [Nq];
		wq_array = new double [Nq];
		gauss(q_array, wq_array, Nq); rangeChange_0_inf(q_array, wq_array, 1000., Nq);
	}

	/* Loops for testing */
	//double P, PP, Q, T;
	////for (int idx_q=0; idx_q<Nq; idx_q++){
	//	double E = Q*Q / MN;
    //  	for (int idx_p=0; idx_p<Np; idx_p++){
	//		P = p_array[idx_p];
	//		for (int idx_pp=0; idx_pp<Np; idx_pp++){
	//			PP = p_array[idx_pp];
	//			
	//		}
	//	}
	//}
	//return 0;

	if (use_premade_antisymmetric_states == true){
		cout << "Reading P123 from file and calculating A123" << endl;
		calculate_antisymmetrization_operator(Np, Nq, Nalpha, &A123, q_array, p_array);

		cout << "Calculating psi_asym from psi_symm using A123" << endl;
		antisymmetrize_state(state_3N_symm_array, state_3N_asym_array, A123, Np, Nq, Nalpha);
	}
	else{
		cout << "Read P123 dimensions from h5-file" << endl;
		int Nalpha_P123 = 0;
		int Np_P123 = 0;
		int Nq_P123 = 0;
		get_h5_P123_dimensions(Nalpha_P123, Np_P123, Nq_P123);
    	P123_array = new double [Np_P123 * Nq_P123 * Nalpha_P123 * Np_P123 * Nq_P123 * Nalpha_P123];

		cout << "Reading P123 from file" << endl;
		read_P123_h5_data_file(P123_array, Nq, q_array, Np, p_array);

		cout << "Starting Faddeev Iterator" << endl;
		calculate_faddeev_convergence(state_3N_asym_array,
									  P123_array,
                     				  Np, p_array, wp_array,
                     				  Nq, q_array, wq_array,
                     				  Nalpha, L_2N, S_2N, J_2N, T_2N, l_3N, two_j_3N,
                     				  two_T_3N, two_J_3N, parity_3N, J_2N_max,
                     				  pot_ptr_nn,
                     				  pot_ptr_np);
	}
	
	/* This function will calculate the 3N c.m. kinetic energy T */
	cout << "Calculating kinetic energy" << endl;
	double kinetic_energy = calculate_3N_kinetic_energy(state_3N_asym_array, state_3N_asym_array, Np, p_array, wp_array, Nq, q_array, wq_array, Nalpha, L_2N, S_2N, J_2N, T_2N, l_3N, two_j_3N);
	
	cout << "Calculating potential energy" << endl;
	double potential_energy = calculate_3N_potential_energy(state_3N_asym_array, state_3N_asym_array, Np, p_array, wp_array, Nq, q_array, wq_array, Nalpha, L_2N, S_2N, J_2N, T_2N, l_3N, two_j_3N, pot_ptr_np, pot_ptr_nn);

	cout << "T: " << kinetic_energy << endl;
	cout << "V: " << potential_energy << endl;
	cout << "H: " << kinetic_energy + potential_energy << endl;

	/* End main body of code here */

	auto program_end= chrono::system_clock::now();
	

	chrono::duration<double> total_time = program_end - program_start;
	cout << "Total run-time:                         " << total_time.count() << endl;

	std::cout << "END OF RUN" << std::endl;
	
    return 0;
}


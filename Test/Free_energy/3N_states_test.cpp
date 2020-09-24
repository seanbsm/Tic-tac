
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

/* Time-keeping modules */
#include <chrono>
#include <ctime>

#include "kinetics.h"
#include "permutation_operators.h"
#include "state_antisymmetrization.h"
#include "Interactions/potential_model.h"
#include "Kai_code_and_data_package/read_psi_Sean.h"

using namespace std;

int main(int argc, char* argv[]){

	auto program_start = chrono::system_clock::now();
	
	/* Start main body of code here */

	/* This function should return a state space basis (for a 2N pair and a free particle) */
	//get_all_states(true);

	/* Quadrature 3N momenta */
	double* p_array  = NULL;
	double* q_array  = NULL;
	double* wp_array = NULL;
	double* wq_array = NULL;

	/* Momentum-representation of 3N state at quadrature nodes in p_array and q_array */
	double* state_3N_symm_array = NULL;
	double* state_3N_asym_array = NULL;

	/* Anti-symmetrization operator */
	double* A123  = NULL;

	/* Quantum numbers of partial-wave expansion in state_3N_array */
    int* L_2N     = NULL;     // pair angular momentum
    int* S_2N     = NULL;     // pair total spin
    int* J_2N     = NULL;     // pair total angular momentum
    int* T_2N     = NULL;     // pair total isospin
    int* l_3N     = NULL;     // three-nucleon angular momentum (?)
    int* two_j_3N = NULL; 	  // three-nucleon total angular momentum x2 (?)

	int Np, Nq, Nalpha;

	potential_model* pot_ptr_np = potential_model::fetch_potential_ptr("Idaho_N3LO", "np");
	potential_model* pot_ptr_nn = potential_model::fetch_potential_ptr("Idaho_N3LO", "nn");

	cout << "Reading states from file" << endl;
	get_all_states(&state_3N_symm_array, &state_3N_asym_array, Np, &p_array, &wp_array, Nq, &q_array, &wq_array, Nalpha, &L_2N, &S_2N, &J_2N, &T_2N, &l_3N, &two_j_3N);

	cout << "Reading P123 from file and calculating A123" << endl;
	calculate_antisymmetrization_operator(Np, Nq, Nalpha, &A123);
	
	cout << "Calculating psi_asym from psi_symm using A123" << endl;
	antisymmetrize_state(state_3N_symm_array, state_3N_asym_array, A123, Np, Nq, Nalpha);
	
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


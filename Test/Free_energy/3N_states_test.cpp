
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

/* Time-keeping modules */
#include <chrono>
#include <ctime>

#include "kinetics.h"
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
	double* state_3N_array = NULL;

	/* Quantum numbers of partial-wave expansion in state_3N_array */
    int* L_2N     = NULL;     // pair angular momentum
    int* S_2N     = NULL;     // pair total spin
    int* J_2N     = NULL;     // pair total angular momentum
    int* T_2N     = NULL;     // pair total isospin
    int* l_3N     = NULL;     // three-nucleon angular momentum (?)
    int* two_j_3N = NULL; // three-nucleon total angular momentum x2 (?)

	int Np, Nq, Nalpha;

	get_all_states(true, &state_3N_array, Np, &p_array, &wp_array, Nq, &q_array, &wq_array, Nalpha, &L_2N, &S_2N, &J_2N, &T_2N, &l_3N, &two_j_3N);

	/* Inner-product test */
    double inner_product = 0;
    for (int idx_p=0; idx_p<Np; idx_p++){
        for (int idx_q=0; idx_q<Nq; idx_q++){
            for (int idx_alpha=0; idx_alpha<Nalpha; idx_alpha++){

                inner_product +=   p_array[idx_p] * p_array[idx_p] * wp_array[idx_p]
                                 * q_array[idx_q] * q_array[idx_q] * wq_array[idx_q]
                                 * state_3N_array[idx_p*Nq*Nalpha + idx_q * Nalpha + idx_alpha]
                                 * state_3N_array[idx_p*Nq*Nalpha + idx_q * Nalpha + idx_alpha];
            }
        }
    }

	cout << inner_product << endl;

	/* This function will construct/read p and q momenta */
	//get_p_and_q_momenta(p_array, q_array, n);

	/* This function will calculate the 3N c.m. kinetic energy T */
	//calculate_3N_kinetic_energy(p_array, q_array, T_array, n);

	/* End main body of code here */

	auto program_end= chrono::system_clock::now();
	

	chrono::duration<double> total_time = program_end - program_start;
	cout << "Total run-time:                         " << total_time.count() << endl;

	std::cout << "END OF RUN" << std::endl;
	
    return 0;
}


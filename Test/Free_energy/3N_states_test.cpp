
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

int main(
	int argc,
	char* argv[]){

	auto program_start = chrono::system_clock::now();
	
	/* Start main body of code here */

	/* This function should return a state space basis (for a 2N pair and a free particle) */
	get_all_states(true);

	const int n = 100;
	double* p_array = new double [n];
	double* q_array = new double [n];
	double* T_array = new double [n];

	/* This function will construct/read p and q momenta */
	get_p_and_q_momenta(p_array, q_array, n);

	/* This function will calculate the 3N c.m. kinetic energy T */
	calculate_3N_kinetic_energy(p_array, q_array, T_array, n);

	/* End main body of code here */

	auto program_end= chrono::system_clock::now();
	

	chrono::duration<double> total_time = program_end - program_start;
	cout << "Total run-time:                         " << total_time.count() << endl;

	std::cout << "END OF RUN" << std::endl;
	
    return 0;
}


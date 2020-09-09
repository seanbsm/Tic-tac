

#include "chiral_N2LOopt.h"

chiral_N2LOopt::chiral_N2LOopt(){
	setup_label_vector_combinations();
	
	chp_set_N2LOopt(&parameters[0]);
}

void chiral_N2LOopt::setup_label_vector_combinations(){
	
	/* 6 parameters are explicitely Tz-dependent */
	int num_lookups = 3 + 3*(labels.size()-3);
	
	state_list.resize( 5*num_lookups );
	state_dep_labels.resize( 2*num_lookups );
	state_dep_idx_lookup.resize( num_lookups );
	
	for (unsigned int k=0; k<3; k++){
		
		/* Append quantum numbers of state */
		state_list[5*k]   = state_list_quantum_numbers[5*k];
		state_list[5*k+1] = state_list_quantum_numbers[5*k + 1];
		state_list[5*k+2] = state_list_quantum_numbers[5*k + 2];
		state_list[5*k+3] = state_list_quantum_numbers[5*k + 3];
		state_list[5*k+4] = state_list_quantum_numbers[5*k + 4];
		
		/* Append index of parameter */
		state_dep_idx_lookup[k] = k;
		
		/* Append label and system for parameter */
		std::string system = system_list[state_list_quantum_numbers[5*k + 4] + 1];
		state_dep_labels[2*k] = labels[k];
		state_dep_labels[2*k+1] = system;
	}
		
	for (unsigned int k=0; k<labels.size()-3; k++){
		for (int Tz=-1; Tz<2; Tz++){
			state_list[5*(3 + 3*k) + 5*(Tz+1) ]    = state_list_quantum_numbers[5*(k+3)];
			state_list[5*(3 + 3*k) + 5*(Tz+1) + 1] = state_list_quantum_numbers[5*(k+3) + 1];
			state_list[5*(3 + 3*k) + 5*(Tz+1) + 2] = state_list_quantum_numbers[5*(k+3) + 2];
			state_list[5*(3 + 3*k) + 5*(Tz+1) + 3] = state_list_quantum_numbers[5*(k+3) + 3];
			state_list[5*(3 + 3*k) + 5*(Tz+1) + 4] = Tz;
			
			/* Append index of parameter
			 * All combinations of Tz point to the same parameter */
			state_dep_idx_lookup[3+3*k + Tz + 1] = k+3;
			
			/* Append label and system for parameter */
			std::string system = system_list[Tz + 1];
			
			state_dep_labels[2*3 + 2*3*k +2*(Tz+1) ]    = labels[k+3];
			state_dep_labels[2*3 + 2*3*k +2*(Tz+1) + 1] = system;
		}
	}
}

void chiral_N2LOopt::setSystem(int Tz_in){
	Tz = -Tz_in;
}

void chiral_N2LOopt::update_parameters(){
	chp_set_N2LOopt(&parameters[0]);
}

void chiral_N2LOopt::setParameter(int index, floatType value){
	parameters[index] = value;
}

void chiral_N2LOopt::retrieveSymmetric(bool &symmetric_in){
	symmetric_in = symmetric;
}

void chiral_N2LOopt::retrieve_state_dep_labels(std::vector<std::string> &state_dep_labels_in){
	state_dep_labels_in = state_dep_labels;
}

void chiral_N2LOopt::retrieve_state_dep_indices(std::vector<int> &state_dep_idx_lookup_in){
	state_dep_idx_lookup_in = state_dep_idx_lookup;
}

void chiral_N2LOopt::retrieve_state_indep_labels(std::vector<std::string> &state_indep_labels_in){
	state_indep_labels_in = state_indep_labels;
}

void chiral_N2LOopt::retrieve_state_indep_indices(std::vector<int> &state_indep_idx_lookup_in){
	state_indep_idx_lookup_in = state_indep_idx_lookup;
}

void chiral_N2LOopt::retrieve_state_list(std::vector<int> &state_list_in){
	state_list_in = state_list;
}

void chiral_N2LOopt::V(double &qi, double &qo, bool coupled, int &S, int &J, int &T, double *Varray){
	
	/* Empty out last set of elements */
	tempVarray[0] = 0;
	tempVarray[1] = 0;
	tempVarray[2] = 0;
	tempVarray[3] = 0;
	tempVarray[4] = 0;
	tempVarray[5] = 0;
	
	int coupledState = coupled;
	
	__idaho_chiral_potential_MOD_chp(&qi, &qo, &coupledState, &S, &J, &T, &Tz, tempVarray);

	/* Change order of array to fit rest of code */
	Varray[0] = tempVarray[0];	//S=0
	Varray[1] = tempVarray[1];	//S=1
	Varray[2] = tempVarray[3];	// Li==Lo<J
	Varray[3] = tempVarray[5];	// Li<Lo
	Varray[4] = tempVarray[4];	// Li>Lo
	Varray[5] = tempVarray[2];	// Li==Lo>J
}

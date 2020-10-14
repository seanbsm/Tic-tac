
#include "chiral_LO_internal.h"

chiral_LO_internal::chiral_LO_internal(){
}

void chiral_LO_internal::setNucleonMass(double M_in){
	M = M_in;
}

void chiral_LO_internal::setPionExchangeClass(int Jmin, int Jmax){
	pionExchange = new OPEP(Jmin, Jmax);
}

void chiral_LO_internal::setSystem(int Tz_in){
	Tz = Tz_in;
}

void chiral_LO_internal::setParameter(int index, floatType value){
	parameters[index] = value;
}

void chiral_LO_internal::update_parameters(){
}

void chiral_LO_internal::retrieveSymmetric(bool &symmetric_in){
	symmetric_in = symmetric;
}

void chiral_LO_internal::retrieve_state_dep_labels(std::vector<std::string> &state_dep_labels_in){
	state_dep_labels_in = state_dep_labels;
}

void chiral_LO_internal::retrieve_state_dep_indices(std::vector<int> &state_dep_idx_lookup_in){
	state_dep_idx_lookup_in = state_dep_idx_lookup;
}

void chiral_LO_internal::retrieve_state_indep_labels(std::vector<std::string> &state_indep_labels_in){
	state_indep_labels_in = state_indep_labels;
}

void chiral_LO_internal::retrieve_state_indep_indices(std::vector<int> &state_indep_idx_lookup_in){
	state_indep_idx_lookup_in = state_indep_idx_lookup;
}

void chiral_LO_internal::retrieve_state_list(std::vector<int> &state_list_in){
	state_list_in = state_list;
}

void chiral_LO_internal::V(double &qi, double &qo, bool coupled, int &S, int &J, int &T, double *Varray){
	
	Varray[0] = 0;
	Varray[1] = 0;
	Varray[2] = 0;
	Varray[3] = 0;
	Varray[4] = 0;
	Varray[5] = 0;
	
	/* Pion-exhange */
	pionExchange->potential(qi, qo, J, Varray);
	
	/* Contact terms */
	if (J==0){
		Varray[0] += parameters[0];
	}
	else if (J==1){
		Varray[2] += parameters[1];
	}
	
	/* Minimal relativity factors */
	double Epi = sqrt(M*M + qi*qi); 	// relativistic energy of in-going particle
	double Epo = sqrt(M*M + qo*qo); 	// relativistic energy of out-going particle
	double relFactor_i = sqrt(M/Epi);	// relativistic factor of in-going particle
	double relFactor_o = sqrt(M/Epo);	// relativistic factor of out-going particle
	
	/* Regulator-functions and Fourier transform constants */
	double temp1 = qi/Lambda;
	double temp2 = temp1*temp1*temp1*temp1*temp1*temp1;	// i.e. we have (qi/Lambda)^6
	double f1 	 = exp(-temp2);
	temp1 = qo/Lambda;
	temp2 = temp1*temp1*temp1*temp1*temp1*temp1;		// i.e. we have (qo/Lambda)^6
	double f2 	 = exp(-temp2);
	double coeff = f1*f2*relFactor_i*relFactor_o/(8*pi*pi*pi);
	
	Varray[0] *= coeff;
	Varray[1] *= coeff;
	Varray[2] *= coeff;
	Varray[3] *= coeff;
	Varray[4] *= coeff;
	Varray[5] *= coeff;
}

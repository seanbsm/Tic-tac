
#ifndef CHIRAL_N2LOOPT_H
#define CHIRAL_N2LOOPT_H

#include <iostream>
#include <string>
#include <vector>

#include "potential_model.h"
#include "../constants.h"

/* External C-function to call fortran script */
extern "C" {
	void chp_set_N2LOopt(double *parameter_array);
    void __idaho_chiral_potential_MOD_chp(double *qi,
			  double *qo,
			  int *coup,
			  int *S,
			  int *J,
			  int *T,
			  int *Tz,
			  double *pot);
}

class chiral_N2LOopt : public potential_model
{
private:
	int Tz;
	
	/* This function generates all Tz alternatives for each
	 * coupling constant. Since there are so many constants,
	 * this function is an easier alternative to writing all
	 * combinations out by hand in the vector state_list */
	void setup_label_vector_combinations();
	
	std::vector<std::string> system_list = {"nn", "np", "pp"};
	
	std::vector<std::string> labels = {"Ct_1S0pp",
									   "Ct_1S0np",
									   "Ct_1S0nn",
									   "Ct_3S1",
									   "C_1S0",
									   "C_3P0",
									   "C_1P1",
									   "C_3P1",
									   "C_3S1",
									   "C_3S1_3D1",
									   "C_3P2"};
									 
	/* Syntax: [S, L1, L2, J, Tz (if applicable)] 
	 * where L1 is the ket ang.mom. and L2 is
	 * the bra ang.mom. */
	const std::vector<int> state_list_quantum_numbers = {0,0,0,0, 1,	// Ct_1S0pp
														 0,0,0,0, 0,	// Ct_1S0np
														 0,0,0,0,-1,	// Ct_1S0nn
														 1,0,0,1, 0,	// Ct_3S1
														 0,0,0,0, 0,	// C_1S0
														 1,1,1,0, 0,	// C_3P0
														 0,1,1,1, 0,	// C_1P1
														 1,1,1,1, 0,	// C_3P1
														 1,0,0,1, 0,	// C_3S1									   
														 1,0,2,1, 0,	// C_3S1_3D1										   
														 1,1,1,2, 0};	// C_3P2
	
	std::vector<std::string> state_dep_labels;
	std::vector<std::string> state_indep_labels = {"gA", 	 	    "np",
												   "c1",		 	"np",
												   "c3",		 	"np",
												   "c4",		 	"np"};
	std::vector<int> state_list;
	std::vector<int> state_dep_idx_lookup;
	std::vector<int> state_indep_idx_lookup = {11,12,13,14};
	
	std::vector<double> parameters = {Ct_1S0pp_n2lo_opt,
									  Ct_1S0np_n2lo_opt,
									  Ct_1S0nn_n2lo_opt,
									  Ct_3S1_n2lo_opt,
									  C_1S0_n2lo_opt,
									  C_3P0_n2lo_opt,
									  C_1P1_n2lo_opt,
									  C_3P1_n2lo_opt,
									  C_3S1_n2lo_opt,
									  C_3S1_3D1_n2lo_opt,
									  C_3P2_n2lo_opt,
									  gA_n2lo_opt,
									  c1_n2lo_opt,
									  c3_n2lo_opt,
									  c4_n2lo_opt};
	
	
	const std::vector<int> indepLabelLookUp = {};
	
	const int numParameters = parameters.size();
	const int numStateIndepParams=0;
	const int numStateDepParams=0;
	const bool symmetric=true;
	
	double *tempVarray = new double [6];
public:
	chiral_N2LOopt();
	
	void retrieveSymmetric(bool &symmetric_in);
	
	void retrieve_state_dep_labels(std::vector<std::string> &state_dep_labels_in);
	void retrieve_state_dep_indices(std::vector<int> &state_dep_idx_lookup_in);
	void retrieve_state_indep_labels(std::vector<std::string> &state_indep_labels_in);
	void retrieve_state_indep_indices(std::vector<int> &state_indep_idx_lookup_in);
	void retrieve_state_list(std::vector<int> &state_list_in);
	
	void setSystem(int Tz_in);
	
	void update_parameters();
	void setParameter(int index, floatType value);
	
	void V(double &qi, double &qo, bool coupled, int &S, int &J, int &T, double *Varray);
};

#endif // CHIRAL_N2LOOPT_H

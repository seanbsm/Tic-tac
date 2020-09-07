
#ifndef CHIRAL_IDAHO_N3LO_H
#define CHIRAL_IDAHO_N3LO_H

#include <iostream>
#include <string>
#include <vector>

#include "potential_model.h"
#include "../constants.h"

/* External C-function to call fortran script */
extern "C" {
	void chp_set_Idaho_N3LO(double *parameter_array);
    void __idaho_chiral_potential_MOD_chp(double *qi,
			  double *qo,
			  int *coup,
			  int *S,
			  int *J,
			  int *T,
			  int *Tz,
			  double *pot);
}

class chiral_Idaho_N3LO : public potential_model
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
									   "C_3P2",
									   "Dh_1S0",
									   "D_1S0",
									   "D_3P0",
									   "D_1P1",
									   "D_3P1",
									   "Dh_3S1",
									   "D_3S1",
									   "D_3D1",
									   "Dh_3S1_3D1",
									   "D_3S1_3D1",
									   "D_1D2",
									   "D_3D2",
									   "D_3P2",
									   "D_3P2_3F2",
									   "D_3D3"};
									 
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
														 1,1,1,2, 0,	// C_3P2										   
														 0,0,0,0, 0,	// Dh_1S0										   
														 0,0,0,0, 0,	// D_1S0										   
														 1,1,1,0, 0,	// D_3P0										   
														 0,1,1,1, 0,	// D_1P1										   
														 1,1,1,1, 0,	// D_3P1										   
														 1,0,0,1, 0,	// Dh_3S1										   
														 1,0,0,1, 0,	// D_3S1										   
														 1,2,2,1, 0,	// D_3D1								   
														 1,0,2,1, 0,	// Dh_3S1_3D1										   
														 1,0,2,1, 0,	// D_3S1_3D1										   
														 0,2,2,2, 0,	// D_1D2										   
														 1,2,2,2, 0,	// D_3D2										   
														 1,1,1,2, 0,	// D_3P2									   
														 1,1,3,2, 0,	// D_3P2_3F2									   
														 1,2,2,3, 0};	// D_3D3
	
	std::vector<std::string> state_dep_labels;
	std::vector<std::string> state_indep_labels = {"gA", 	 	    "np",
												   "c1",		 	"np",
												   "c2",		 	"np",
												   "c3",		 	"np",
												   "c4",		 	"np",
												   "d1_plus_d2", 	"np",
												   "d3",		 	"np",
												   "d5",		 	"np",
												   "d14_minus_d15", "np"};
	std::vector<int> state_list;
	std::vector<int> state_dep_idx_lookup;
	std::vector<int> state_indep_idx_lookup = {26,27,28,29,30,31,32,33,34};
	
	std::vector<double> parameters = {Ct_1S0pp_idaho_n3lo,
									  Ct_1S0np_idaho_n3lo,
									  Ct_1S0nn_idaho_n3lo,
									  Ct_3S1_idaho_n3lo,
									  C_1S0_idaho_n3lo,
									  C_3P0_idaho_n3lo,
									  C_1P1_idaho_n3lo,
									  C_3P1_idaho_n3lo,
									  C_3S1_idaho_n3lo,
									  C_3S1_3D1_idaho_n3lo,
									  C_3P2_idaho_n3lo,
									  Dh_1S0_idaho_n3lo,
									  D_1S0_idaho_n3lo,
									  D_3P0_idaho_n3lo,
									  D_1P1_idaho_n3lo,
									  D_3P1_idaho_n3lo,
									  Dh_3S1_idaho_n3lo,
									  D_3S1_idaho_n3lo,
									  D_3D1_idaho_n3lo,
									  Dh_3S1_3D1_idaho_n3lo,
									  D_3S1_3D1_idaho_n3lo,
									  D_1D2_idaho_n3lo,
									  D_3D2_idaho_n3lo,
									  D_3P2_idaho_n3lo,
									  D_3P2_3F2_idaho_n3lo,
									  D_3D3_idaho_n3lo,
									  gA_idaho_n3lo,
									  c1_idaho_n3lo,
									  c2_idaho_n3lo,
									  c3_idaho_n3lo,
									  c4_idaho_n3lo,
									  d1_plus_d2_idaho_n3lo,
									  d3_idaho_n3lo,
									  d5_idaho_n3lo,
									  d14_minus_d15_idaho_n3lo};
	
	
	const std::vector<int> indepLabelLookUp = {};
	
	const int numParameters = parameters.size();
	const int numStateIndepParams=0;
	const int numStateDepParams=0;
	const bool symmetric=true;
	
	double *tempVarray = new double [6];
public:
	chiral_Idaho_N3LO();
	
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

#endif // CHIRAL_IDAHO_N3LO_H

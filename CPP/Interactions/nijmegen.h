
#ifndef NIJMEGEN_H
#define NIJMEGEN_H

#include <iostream>
#include <string>
#include <vector>

#include "potential_model.h"
#include "../constants.h"

/* External C-function to call fortran script */
extern "C" {
    void nijmegen_fort_interface(double *qi,
			  double *qo,
			  int *coup,
			  int *S,
			  int *J,
			  int *T,
			  int *Tz,
			  double *pot);
}

class nijmegen : public potential_model
{
private:
	double M;
	int Tz;
	
	std::vector<floatType> parameters = {};
	
	const std::vector<std::string> state_dep_labels = {};
	const std::vector<std::string> state_indep_labels;
	
	const std::vector<int> state_dep_idx_lookup = {};
	const std::vector<int> state_indep_idx_lookup;
	
	/* Syntax: [S, L1, L2, J, Tz, paramIndex] 
	 * where L1 is the ket ang.mom. and L2 is
	 * the bra ang.mom. */
	const std::vector<int> state_list  = {};
	
	const std::vector<int> indepLabelLookUp = {};
	
	const bool symmetric = true;

	double *tempVarray = new double [6];
public:
	nijmegen();
	
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

#endif // NIJMEGEN_H

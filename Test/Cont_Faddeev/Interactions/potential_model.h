
#ifndef POTENTIAL_MODEL_H
#define POTENTIAL_MODEL_H

#include <iostream>
#include <string>
#include <vector>

#include "../type_defs.h"

class potential_model
{
public:
	potential_model();
	static potential_model *fetch_potential_ptr(std::string model, std::string system);
	
	virtual void retrieveSymmetric(bool &symmetric_in) 		  = 0;
	
	virtual void retrieve_state_dep_labels(std::vector<std::string> &state_dep_labels_in)	  = 0;
	virtual void retrieve_state_dep_indices(std::vector<int> &state_dep_idx_lookup_in)   	  = 0;
	virtual void retrieve_state_indep_labels(std::vector<std::string> &state_indep_labels_in) = 0;
	virtual void retrieve_state_indep_indices(std::vector<int> &state_indep_idx_lookup_in)    = 0;
	virtual void retrieve_state_list(std::vector<int> &state_list_in) 						  = 0;
	
	virtual void setSystem(int Tz) 						   										 = 0;
	virtual void setParameter(int index, floatType value)  										 = 0;
	virtual void update_parameters()  															 = 0;
	virtual void V(double &qi, double &qo, bool coupled, int &S, int &J, int &T, double *Varray) = 0;
};

#endif // POTENTIAL_MODEL_H


#ifndef CHIRAL_LO_INTERNAL_H
#define CHIRAL_LO_INTERNAL_H

#include <iostream>
#include <string>
#include <vector>

#include "potential_model.h"
#include "OPEP.h"
#include "../constants.h"

class chiral_LO_internal : public potential_model
{
private:
	double M;
	int Tz;
	
	std::vector<floatType> parameters = {C1S0,
										 C3S1};
	
	const std::vector<std::string> state_dep_labels = {"C1S0", "nn",
													   "C1S0", "np",
													   "C1S0", "pp",
													   "C3S1", "nn",
													   "C3S1", "np",
													   "C3S1", "pp"};
	const std::vector<std::string> state_indep_labels;
	
	const std::vector<int> state_dep_idx_lookup = {0,0,0,
												   1,1,1};
	const std::vector<int> state_indep_idx_lookup;
	
	/* Syntax: [S, L1, L2, J, Tz, paramIndex] 
	 * where L1 is the ket ang.mom. and L2 is
	 * the bra ang.mom. */
	const std::vector<int> state_list  = {0,0,0,0, 1,
										  0,0,0,0, 0,
										  0,0,0,0,-1,
										  1,0,0,1, 1,
										  1,0,0,1, 0,
										  1,0,0,1,-1};
	
	const std::vector<int> indepLabelLookUp = {};
	
	const bool symmetric = true;
	
	OPEP *pionExchange;
public:
	chiral_LO_internal();
	
	void retrieveSymmetric(bool &symmetric_in);
	
	void retrieve_state_dep_labels(std::vector<std::string> &state_dep_labels_in);
	void retrieve_state_dep_indices(std::vector<int> &state_dep_idx_lookup_in);
	void retrieve_state_indep_labels(std::vector<std::string> &state_indep_labels_in);
	void retrieve_state_indep_indices(std::vector<int> &state_indep_idx_lookup_in);
	void retrieve_state_list(std::vector<int> &state_list_in);
	
	void setPionExchangeClass(int Jmin, int Jmax);
	void setNucleonMass(double M_in);
	void setSystem(int Tz_in);
	
	void update_parameters();
	void setParameter(int index, floatType value);
	
	void V(double &qi, double &qo, bool coupled, int &S, int &J, int &T, double *Varray);
};

#endif // CHIRAL_LO_INTERNAL_H

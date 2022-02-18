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
	
	std::vector<floatType> parameters_LO = {C1S0,
										    C3S1};
	
	void setPionExchangeClass(int Jmin, int Jmax);
	void setNucleonMass(double M_in);
	
	OPEP *pionExchange;
public:
	chiral_LO_internal(double MN_in, int Jmin, int Jmax);
	
	void update_parameters(double* parameters);
	
	void V(double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray);
};

#endif // CHIRAL_LO_INTERNAL_H


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
	static potential_model *fetch_potential_ptr(run_params run_parameters);
	
	virtual void update_parameters(double* parameters) = 0;
	
	virtual void V(double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray) = 0;
};

#endif // POTENTIAL_MODEL_H


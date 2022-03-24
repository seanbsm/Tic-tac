
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
	
	virtual void first_parameter_sampling(bool statement) = 0;
	virtual void update_parameters(double* parameters) = 0;
	virtual void setup_store_matrices(double* p_mesh, int Np, bool coupled, int &S, int &J, int &T, int &Tz) = 0;
	
	virtual void V(int i, int j, double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray) = 0;
};

#endif // POTENTIAL_MODEL_H


#ifndef MALFLIET_TJON_H
#define MALFLIET_TJON_H

#include <iostream>
#include <string>
#include <vector>

#include "potential_model.h"
#include "../constants.h"

class malfliet_tjon : public potential_model
{
public:
	malfliet_tjon();

	void first_parameter_sampling(bool statement);
	void update_parameters(double* parameters);
	void setup_store_matrices(double* p_mesh, int Np, bool coupled, int &S, int &J, int &T, int &Tz);
	
	void V(int i, int j, double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray);
};

#endif // MALFLIET_TJON_H

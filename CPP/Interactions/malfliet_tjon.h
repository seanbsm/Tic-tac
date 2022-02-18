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

	void update_parameters(double* parameters);
	
	void V(double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray);
};

#endif // MALFLIET_TJON_H

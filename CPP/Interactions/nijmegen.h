
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
public:
	nijmegen();

	void first_parameter_sampling(bool statement);
	void update_parameters(double* parameters);
	void setup_store_matrices(double* p_mesh, int Np, bool coupled, int &S, int &J, int &T, int &Tz);
	
	void V(int i, int j, double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray);
};

#endif // NIJMEGEN_H


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
	double *tempVarray = new double [6];
public:
	nijmegen();

	void update_parameters(double* parameters);
	
	void V(double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray);
};

#endif // NIJMEGEN_H

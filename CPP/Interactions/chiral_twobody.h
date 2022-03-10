
#ifndef CHIRAL_TWOBODY_H
#define CHIRAL_TWOBODY_H

#include <iostream>
#include <string>
#include <vector>

#include "potential_model.h"
#include "../error_management.h"
#include "../constants.h"

/* External C-function to call fortran script */
extern "C" {
	void chp_set_ispot_lo(double *parameter_array);
	void chp_set_ispot_nlo(double *parameter_array);
	void chp_set_ispot_n2lo(double *parameter_array);
	void chp_set_ispot_n3lo(double *parameter_array);
	
    void __idaho_chiral_potential_MOD_chp(double *qi,
			  double *qo,
			  int *coup,
			  int *S,
			  int *J,
			  int *T,
			  int *Tz,
			  double *pot);
}

class chiral_twobody : public potential_model
{
private:
	std::string model;
public:
	chiral_twobody();

	void call_preset(std::string input_model);

	void update_parameters(double* parameters);
	
	void V(double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray);
};

#endif // CHIRAL_TWOBODY_H

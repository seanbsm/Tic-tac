
#ifndef CHIRAL_N2LOOPT_H
#define CHIRAL_N2LOOPT_H

#include <iostream>
#include <string>
#include <vector>

#include "potential_model.h"
#include "../constants.h"

/* External C-function to call fortran script */
extern "C" {
	void chp_set_N2LOopt(double *parameter_array);
    void __idaho_chiral_potential_MOD_chp(double *qi,
			  double *qo,
			  int *coup,
			  int *S,
			  int *J,
			  int *T,
			  int *Tz,
			  double *pot);
}

class chiral_N2LOopt : public potential_model
{
private:
	std::vector<double> parameters_n2lo_opt = {Ct_1S0pp_n2lo_opt,
									  		   Ct_1S0np_n2lo_opt,
									  		   Ct_1S0nn_n2lo_opt,
									  		   Ct_3S1_n2lo_opt,
									  		   C_1S0_n2lo_opt,
									  		   C_3P0_n2lo_opt,
									  		   C_1P1_n2lo_opt,
									  		   C_3P1_n2lo_opt,
									  		   C_3S1_n2lo_opt,
									  		   C_3S1_3D1_n2lo_opt,
									  		   C_3P2_n2lo_opt,
									  		   gA_n2lo_opt,
									  		   c1_n2lo_opt,
									  		   c3_n2lo_opt,
									  		   c4_n2lo_opt};
									  
	double *tempVarray = new double [6];
public:
	chiral_N2LOopt();

	void update_parameters(double* parameters);
	
	void V(double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray);
};

#endif // CHIRAL_N2LOOPT_H

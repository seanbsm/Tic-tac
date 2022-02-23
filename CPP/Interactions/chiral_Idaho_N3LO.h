
#ifndef CHIRAL_IDAHO_N3LO_H
#define CHIRAL_IDAHO_N3LO_H

#include <iostream>
#include <string>
#include <vector>

#include "potential_model.h"
#include "../constants.h"

/* External C-function to call fortran script */
extern "C" {
	void chp_set_Idaho_N3LO(double *parameter_array);
    void __idaho_chiral_potential_MOD_chp(double *qi,
			  double *qo,
			  int *coup,
			  int *S,
			  int *J,
			  int *T,
			  int *Tz,
			  double *pot);
}

class chiral_Idaho_N3LO : public potential_model
{
private:
	
	std::vector<double> parameters_idaho_n3lo = {Ct_1S0pp_idaho_n3lo,
												 Ct_1S0np_idaho_n3lo,
												 Ct_1S0nn_idaho_n3lo,
												 Ct_3S1_idaho_n3lo,
												 C_1S0_idaho_n3lo,
												 C_3P0_idaho_n3lo,
												 C_1P1_idaho_n3lo,
												 C_3P1_idaho_n3lo,
												 C_3S1_idaho_n3lo,
												 C_3S1_3D1_idaho_n3lo,
												 C_3P2_idaho_n3lo,
												 Dh_1S0_idaho_n3lo,
												 D_1S0_idaho_n3lo,
												 D_3P0_idaho_n3lo,
												 D_1P1_idaho_n3lo,
												 D_3P1_idaho_n3lo,
												 Dh_3S1_idaho_n3lo,
												 D_3S1_idaho_n3lo,
												 D_3D1_idaho_n3lo,
												 Dh_3S1_3D1_idaho_n3lo,
												 D_3S1_3D1_idaho_n3lo,
												 D_1D2_idaho_n3lo,
												 D_3D2_idaho_n3lo,
												 D_3P2_idaho_n3lo,
												 D_3P2_3F2_idaho_n3lo,
												 D_3D3_idaho_n3lo,
												 gA_idaho_n3lo,
												 c1_idaho_n3lo,
												 c2_idaho_n3lo,
												 c3_idaho_n3lo,
												 c4_idaho_n3lo,
												 d1_plus_d2_idaho_n3lo,
												 d3_idaho_n3lo,
												 d5_idaho_n3lo,
												 d14_minus_d15_idaho_n3lo};
public:
	chiral_Idaho_N3LO();
	
	void update_parameters(double* parameters);
	
	void V(double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray);
};

#endif // CHIRAL_IDAHO_N3LO_H

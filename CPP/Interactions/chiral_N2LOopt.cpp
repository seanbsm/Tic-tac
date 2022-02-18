

#include "chiral_N2LOopt.h"

chiral_N2LOopt::chiral_N2LOopt(){
	chp_set_N2LOopt(&parameters_n2lo_opt[0]);
}

void chiral_N2LOopt::update_parameters(double* parameters){
	chp_set_N2LOopt(parameters);
}

void chiral_N2LOopt::V(double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray){
	
	/* Empty out last set of elements */
	tempVarray[0] = 0;
	tempVarray[1] = 0;
	tempVarray[2] = 0;
	tempVarray[3] = 0;
	tempVarray[4] = 0;
	tempVarray[5] = 0;
	
	int coupledState = coupled;
    int Tz_reverse = -Tz;
	
	__idaho_chiral_potential_MOD_chp(&qi, &qo, &coupledState, &S, &J, &T, &Tz_reverse, tempVarray);

	/* Change order of array to fit rest of code */
	Varray[0] = tempVarray[0];	//S=0
	Varray[1] = tempVarray[1];	//S=1
	Varray[2] = tempVarray[3];	// Li==Lo<J
	Varray[3] = tempVarray[5];	// Li<Lo
	Varray[4] = tempVarray[4];	// Li>Lo
	Varray[5] = tempVarray[2];	// Li==Lo>J
}

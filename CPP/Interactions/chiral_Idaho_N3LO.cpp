

#include "chiral_Idaho_N3LO.h"

chiral_Idaho_N3LO::chiral_Idaho_N3LO(){
	chp_set_Idaho_N3LO(&parameters_idaho_n3lo[0]);
}

void chiral_Idaho_N3LO::update_parameters(double* parameters){
	chp_set_Idaho_N3LO(parameters);
}

void chiral_Idaho_N3LO::V(double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray){
	
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
	Varray[0] = tempVarray[0];
	Varray[1] = tempVarray[1];
	Varray[2] = tempVarray[3];
	Varray[3] = tempVarray[5];
	Varray[4] = tempVarray[4];
	Varray[5] = tempVarray[2];
}

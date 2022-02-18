

#include "nijmegen.h"

nijmegen::nijmegen(){
}

void nijmegen::update_parameters(double* parameters){}

void nijmegen::V(double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray){
	
	/* Empty out last set of elements */
	tempVarray[0] = 0;
	tempVarray[1] = 0;
	tempVarray[2] = 0;
	tempVarray[3] = 0;
	tempVarray[4] = 0;
	tempVarray[5] = 0;
	
	int coupledState = coupled;
    int Tz_reverse = -Tz;
	
	nijmegen_fort_interface(&qi, &qo, &coupledState, &S, &J, &T, &Tz_reverse, tempVarray);

	/* Change order of array to fit rest of code */
	Varray[0] = tempVarray[0];
	Varray[1] = tempVarray[1];
	Varray[2] = tempVarray[3];
	Varray[3] = tempVarray[5];
	Varray[4] = tempVarray[4];
	Varray[5] = tempVarray[2];
}

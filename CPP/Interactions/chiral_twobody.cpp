
#include "chiral_twobody.h"

chiral_twobody::chiral_twobody(){
}

void chiral_twobody::call_preset(std::string input_model){
	model = input_model;

	double* parameters = NULL;
	if (input_model=="IS_LO"){
		parameters = new double [2];
	}
	else if (input_model=="IS_NLO"){
		parameters = new double [13];
	}
	else if (input_model=="IS_N2LO"){
		parameters = new double [16];
	}
	else if (input_model=="IS_N3LO"){
		parameters = new double [31];
	}
	else{
		raise_error("Requested chiral potential " + input_model + " doesn't exist in program.");
	}

	update_parameters(parameters);

	delete [] parameters;
}

void chiral_twobody::update_parameters(double* parameters){
	if (model=="IS_LO"){
		chp_set_ispot_lo(parameters);
	}
	else if (model=="IS_NLO"){
		chp_set_ispot_nlo(parameters);
	}
	else if (model=="IS_N2LO"){
		chp_set_ispot_n2lo(parameters);
	}
	else if (model=="IS_N3LO"){
		chp_set_ispot_n3lo(parameters);
	}
	else{
		raise_error("Requested chiral potential doesn't exist in program.");
	}
}

void chiral_twobody::V(double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray){
	
	double tempVarray [6];
	
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

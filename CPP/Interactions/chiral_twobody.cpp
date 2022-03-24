
#include "chiral_twobody.h"

chiral_twobody::chiral_twobody(){
}

void chiral_twobody::call_preset(std::string input_model){
	model = input_model;

	double* parameters = NULL;
	if (input_model=="IS_LO"){
		parameters 	    = new double [2];
		parameters_copy = new double [2];
	}
	else if (input_model=="IS_NLO"){
		parameters 	    = new double [13];
		parameters_copy = new double [13];
	}
	else if (input_model=="IS_N2LO"){
		parameters  	= new double [16];
		parameters_copy = new double [16];
	}
	else if (input_model=="IS_N3LO"){
		parameters 	    = new double [31];
		parameters_copy = new double [31];
	}
	else{
		raise_error("Requested chiral potential " + input_model + " doesn't exist in program.");
	}

	update_parameters(parameters);

	delete [] parameters;
}

void chiral_twobody::first_parameter_sampling(bool statement){
	first_sampling = statement;
}

void chiral_twobody::update_parameters(double* parameters){
	if (model=="IS_LO"){
		chp_set_ispot_lo(parameters);
		for(int i=0;i<2;i++){parameters_copy[i]=parameters[i];}	// Copy parameter array
	}
	else if (model=="IS_NLO"){
		chp_set_ispot_nlo(parameters);
		for(int i=0;i<13;i++){parameters_copy[i]=parameters[i];}	// Copy parameter array
	}
	else if (model=="IS_N2LO"){
		chp_set_ispot_n2lo(parameters);
		for(int i=0;i<16;i++){parameters_copy[i]=parameters[i];}	// Copy parameter array
	}
	else if (model=="IS_N3LO"){
		chp_set_ispot_n3lo(parameters);
		for(int i=0;i<31;i++){parameters_copy[i]=parameters[i];}	// Copy parameter array
	}
	else{
		raise_error("Requested chiral potential doesn't exist in program.");
	}
}

void chiral_twobody::set_run_parameters(run_params run_parameters_in){
	run_parameters = run_parameters_in;

	if (V_2L2PE_allocated==false && model=="IS_N3LO"){
		/* Figure out momentum-mesh dimension */
		if (run_parameters.midpoint_approx==true){
			Np = run_parameters.Np_WP;
		}
		else{
			Np = run_parameters.Np_WP * run_parameters.Np_per_WP;
		}
		/* Number of Tz-options */
		NTz = 3; // Tz=-1,0,1
		int tot_size = NTz * Np*Np * (num_2L2PE_params+1) * 6;
		V_2L2PE_array = new double [tot_size];
		for (int i=0; i<tot_size; i++){
			V_2L2PE_array[i] = 0;
		}
		V_2L2PE_allocated = true;
	}
}

void chiral_twobody::setup_store_matrices(double* p_mesh, int Np_in, bool coupled, int &S, int &J, int &T, int &Tz){
	/* Construct interaction for first time */
	if (first_sampling && model=="IS_N3LO"){

		double* parameters = new double [31];
		double* tempVarray = new double [6];
		int coupledState = coupled;
    	int Tz_reverse = -Tz;
		
		/* Set constant term */
		for (int k=0; k<31; k++){
			parameters[k] = 0;
		}
		chp_set_ispot_n3lo(parameters);
		for (int i=0; i<Np; i++){
			double qo = p_mesh[i];
			for (int j=0; j<Np; j++){
				double qi = p_mesh[j];
				__idaho_chiral_potential_MOD_chp(&qi, &qo, &coupledState, &S, &J, &T, &Tz_reverse, tempVarray);
				//tempVarray[0] = 5.1; tempVarray[1] = 5.2; tempVarray[2] = 5.4; tempVarray[3] = 5.6; tempVarray[4] = 5.1; tempVarray[5] = 5.1;	// I used this to omit the calculation time while debugging indexing

				int ij_arr_idx = ((Tz+1)*Np*Np + i*Np + j) * (num_2L2PE_params+1) * 6;
				/* Add constant term */
				for (int k=0; k<6; k++){
					V_2L2PE_array[ij_arr_idx + k] += tempVarray[k];	
					tempVarray[k] = 0;
				}
			}
		}
		
		
		/* Set interaction terms */
		for (int p=1; p<num_2L2PE_params+1; p++){
			parameters[3+p] = 1.0;
			chp_set_ispot_n3lo(parameters);
			
			for (int i=0; i<Np; i++){
				double qo = p_mesh[i];
				for (int j=0; j<Np; j++){
					double qi = p_mesh[j];
					__idaho_chiral_potential_MOD_chp(&qi, &qo, &coupledState, &S, &J, &T, &Tz_reverse, tempVarray);
					//tempVarray[0] = 5.1; tempVarray[1] = 5.2; tempVarray[2] = 5.4; tempVarray[3] = 5.6; tempVarray[4] = 5.1; tempVarray[5] = 5.1;	// I used this to omit the calculation time while debugging indexing

					int ij_arr_idx = ((Tz+1)*Np*Np + i*Np + j) * (num_2L2PE_params+1) * 6;
					/* Add coefficient terms and subtract constant term */
					for (int k=0; k<6; k++){
						if (tempVarray[k]!=0){
							V_2L2PE_array[ij_arr_idx + p*num_2L2PE_params + k] += tempVarray[k] - V_2L2PE_array[ij_arr_idx + k];
						}
						tempVarray[k] = 0;
					}
				}
			}
			parameters[3+p] = 0.0;
		}

		/* Delete temporary array */
		delete [] parameters;
		delete [] tempVarray;

		/* Restore parameters in chp-preset */
		chp_set_ispot_n3lo(parameters_copy);
	}
}

/* This function uses std:map to minimize the number of 2-loop 2-pion calculations.
 * Probably, std:unordered_map (i.e. a hash map) would be faster, but requires a custom hash-function. */
void chiral_twobody::calculate_2loop_2PE(int i, int j, double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *tempVarray){
	/* Since Tz,i,j are the outermost indices of V_2L2PE_array, they define a constant steplength throughout this function */
	int ij_arr_idx = ((Tz+1)*Np*Np + i*Np + j) * (num_2L2PE_params+1) * 6;

	/* Use stored interaction to fill tempVarray */
	for (int k=0; k<6; k++){
		/* Add constant term */
		//tempVarray[k] = V_2L2PE_array[ij_arr_idx + k];	// THIS IS ALREADY ADDED IN chiral_twobody::V(...){...}
		for (int p=1; p<num_2L2PE_params+1; p++){
			/* Add coefficient terms */
			tempVarray[k] += parameters_copy[3+p] * V_2L2PE_array[ij_arr_idx + p*num_2L2PE_params + k] / 100.;	// LECs are given in 10^4 GeV^-2, but we use MeV^-2
		}
	}
}

void chiral_twobody::V(int i, int j, double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray){
	
	double tempVarray [6] = {0};
	
	int coupledState = coupled;
    int Tz_reverse = -Tz;

	if (model=="IS_N3LO" && false){
		/* Set 2-loop 2-pion exchange variables to zero and use instead pre-stored interactions */
		double* parameters = new double [31];
		for(int i=0;i<31;i++){parameters[i]=parameters_copy[i];} // Copy parameter array
		for(int i=4;i<5;i++){parameters[i]=0;} // Set 2-loop 2-pion exchange variables to zero
		chp_set_ispot_n3lo(parameters);
		delete [] parameters;
		/* Evaluate other diagrams */
		__idaho_chiral_potential_MOD_chp(&qi, &qo, &coupledState, &S, &J, &T, &Tz_reverse, tempVarray);
		/* BACKGROUND */
		double tempVarray2 [6] = {0};
		parameters = new double [31];
		for(int i=0;i<31;i++){parameters[i]=0;} // Copy parameter array
		chp_set_ispot_n3lo(parameters);
		delete [] parameters;
		__idaho_chiral_potential_MOD_chp(&qi, &qo, &coupledState, &S, &J, &T, &Tz_reverse, tempVarray2);
		/* Use pre-stored interatctions for 2-loop 2-pion exchange */
		parameters = new double [31];
		for(int i=0;i<31;i++){parameters[i]=0;} // Copy parameter array
		for(int i=4;i<5;i++){
			double tempVarray1 [6] = {0};
			
			parameters[i]=1.;//parameters_copy[i]; // Set 2-loop 2-pion exchange variables to zero
			chp_set_ispot_n3lo(parameters);
			__idaho_chiral_potential_MOD_chp(&qi, &qo, &coupledState, &S, &J, &T, &Tz_reverse, tempVarray1);
			for (int j=0; j<6;j++){
				tempVarray[j] += parameters_copy[i]*tempVarray1[j]-tempVarray2[j];
				tempVarray1[j] = 0;
			}
			parameters[i] = 0;
		}
		delete [] parameters;
		
	}
	else{
		__idaho_chiral_potential_MOD_chp(&qi, &qo, &coupledState, &S, &J, &T, &Tz_reverse, tempVarray);
	}

	/* Change order of array to fit rest of code */
	Varray[0] = tempVarray[0];
	Varray[1] = tempVarray[1];
	Varray[2] = tempVarray[3];
	Varray[3] = tempVarray[5];
	Varray[4] = tempVarray[4];
	Varray[5] = tempVarray[2];
}

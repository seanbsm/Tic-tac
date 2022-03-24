
#ifndef CHIRAL_TWOBODY_H
#define CHIRAL_TWOBODY_H

#include <iostream>
#include <string>
#include <vector>
#include <map>		// Used in the case of N3LO

#include "potential_model.h"
#include "../type_defs.h"
#include "../make_pw_symm_states.h"
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
	
	/* This array copies the input of the latest call to update_parameters(). */
	double* parameters_copy;

	int Np;
	int NTz;

	bool first_sampling;

	bool    V_2L2PE_allocated = false;
	double* V_2L2PE_array = NULL;

	run_params run_parameters;
	const int num_2L2PE_params = 4;	// 2L2PE = "two-loop two-pion exchange"

	void calculate_2loop_2PE(int i, int j, double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *tempVarray);
public:
	chiral_twobody();

	void call_preset(std::string input_model);
	void set_run_parameters(run_params run_parameters_in);
	void setup_store_matrices(double* p_mesh, int Np, bool coupled, int &S, int &J, int &T, int &Tz);

	void first_parameter_sampling(bool statement);
	void update_parameters(double* parameters);
	
	void V(int i, int j, double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray);
};

#endif // CHIRAL_TWOBODY_H

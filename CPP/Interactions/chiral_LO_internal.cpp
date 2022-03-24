
#include "chiral_LO_internal.h"

chiral_LO_internal::chiral_LO_internal(double MN_in, int Jmin, int Jmax){
	setPionExchangeClass(Jmin, Jmax);
	setNucleonMass(MN_in);
}

void chiral_LO_internal::setNucleonMass(double M_in){
	M = M_in;
}

void chiral_LO_internal::setPionExchangeClass(int Jmin, int Jmax){
	pionExchange = new OPEP(Jmin, Jmax);
}

void chiral_LO_internal::first_parameter_sampling(bool statement){
}

void chiral_LO_internal::update_parameters(double* parameters){
}

void chiral_LO_internal::setup_store_matrices(double* p_mesh, int Np, bool coupled, int &S, int &J, int &T, int &Tz){
}

void chiral_LO_internal::V(int i, int j, double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray){
	
	Varray[0] = 0;
	Varray[1] = 0;
	Varray[2] = 0;
	Varray[3] = 0;
	Varray[4] = 0;
	Varray[5] = 0;
	
	/* Pion-exhange */
	pionExchange->potential(qi, qo, J, Varray);
	
	/* Contact terms */
	if (J==0){
		Varray[0] += parameters_LO[0];
	}
	else if (J==1){
		Varray[2] += parameters_LO[1];
	}
	
	/* Minimal relativity factors */
	double Epi = sqrt(M*M + qi*qi); 	// relativistic energy of in-going particle
	double Epo = sqrt(M*M + qo*qo); 	// relativistic energy of out-going particle
	double relFactor_i = sqrt(M/Epi);	// relativistic factor of in-going particle
	double relFactor_o = sqrt(M/Epo);	// relativistic factor of out-going particle
	
	/* Regulator-functions and Fourier transform constants */
	double temp1 = qi/Lambda;
	double temp2 = temp1*temp1*temp1*temp1*temp1*temp1;	// i.e. we have (qi/Lambda)^6
	double f1 	 = exp(-temp2);
	temp1 = qo/Lambda;
	temp2 = temp1*temp1*temp1*temp1*temp1*temp1;		// i.e. we have (qo/Lambda)^6
	double f2 	 = exp(-temp2);
	double coeff = f1*f2*relFactor_i*relFactor_o/(8*pi*pi*pi);
	
	Varray[0] *= coeff;
	Varray[1] *= coeff;
	Varray[2] *= coeff;
	Varray[3] *= coeff;
	Varray[4] *= coeff;
	Varray[5] *= coeff;
}

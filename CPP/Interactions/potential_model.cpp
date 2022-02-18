
#include "potential_model.h"

#include "../constants.h"

/* Include present potentials */
#include "chiral_LO_internal.h"
#include "chiral_N2LOopt.h"
#include "chiral_Idaho_N3LO.h"
#include "malfliet_tjon.h"
#include "nijmegen.h"


potential_model::potential_model(){
}

potential_model *potential_model::fetch_potential_ptr(std::string model){
	
	if (model=="LO_internal"){
		chiral_LO_internal *pot_ptr = new chiral_LO_internal(MN, 0, 100);
		
		std::cout << "Beware that the chiral_LO_internal potential precalculates the Legendre polynomials for a given range [Jmin, Jmax] of J-values. \n"
				  << "This range is by default set to Jmin=0 and Jmax=100. If you really require calculations for Jmax>100, you can change this \n"
				  << "BEFORE compilation in potential_model.cpp " << std::endl;

		return pot_ptr;
	}
	else if (model=="N2LOopt"){
		chiral_N2LOopt *pot_ptr = new chiral_N2LOopt();
		return pot_ptr;
	}
	else if (model=="Idaho_N3LO"){
		chiral_Idaho_N3LO *pot_ptr = new chiral_Idaho_N3LO();
		return pot_ptr;
	}
	else if (model=="malfliet_tjon"){
		malfliet_tjon* pot_ptr = new malfliet_tjon();
		return pot_ptr;
	}
	else if (model=="nijmegen"){
		nijmegen* pot_ptr = new nijmegen();
		return pot_ptr;
	}
	//else if (model=="chiral-twobody"){
	//	chp* pot_ptr = new chp();
	//	return pot_ptr;
	//}
	else{
		std::cout << "Invalid potential model entered. Exiting ..." << std::endl;
		exit(-1);
	}
}

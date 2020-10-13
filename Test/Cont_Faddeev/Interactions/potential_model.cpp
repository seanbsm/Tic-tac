
#include "potential_model.h"

#include "../constants.h"

/* Include present potentials */
#include "chiral_LO_internal.h"
#include "chiral_N2LOopt.h"
#include "chiral_Idaho_N3LO.h"


potential_model::potential_model(){
}

int find_Tz_from_system(std::string system){
	int Tz;
	if (system=="np"){
		Tz = 0;
	}
	else if (system=="pp"){
		Tz = 1;
	}
	else if (system=="nn"){
		Tz = -1;
	}
	else{
		std::cout << "Invalid system entered. Exiting ..." << std::endl;
		exit(0);
	}
	
	return Tz;
}

potential_model *potential_model::fetch_potential_ptr(std::string model, std::string system){
	
	int Tz = find_Tz_from_system(system);
	
	if (model=="LO_internal"){
		chiral_LO_internal *pot_ptr = new chiral_LO_internal();
		
		pot_ptr->setSystem(Tz);
		pot_ptr->setNucleonMass(MN);

		std::cout << "Beware that the chiral_LO_internal potential precalculates the Legendre polynomials for a given range [Jmin, Jmax] of J-values. \n"
				  << "This range is by default set to Jmin=0 and Jmax=100. If you really require calculations for Jmax>100, you can change this \n"
				  << "BEFORE compilation in potential_model.cpp " << std::endl;
		pot_ptr->setPionExchangeClass(0, 100);

		return pot_ptr;
	}
	else if (model=="N2LOopt"){
		chiral_N2LOopt *pot_ptr = new chiral_N2LOopt();
		
		pot_ptr->setSystem(Tz);
		
		return pot_ptr;
	}
	else if (model=="Idaho_N3LO"){
		chiral_Idaho_N3LO *pot_ptr = new chiral_Idaho_N3LO();
		
		pot_ptr->setSystem(Tz);
		
		return pot_ptr;
	}
	else{
		std::cout << "Invalid potential model entered. Exiting ..." << std::endl;
		exit(-1);
	}
}

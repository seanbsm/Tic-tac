
#include "potential_model.h"

/* Include present potentials */
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
	
	if (model=="N2LOopt"){
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

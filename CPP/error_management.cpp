
#include "error_management.h"

void raise_error (std::string error_cause){
	if (error_cause == "problem_reading_data"){
		std::cout <<
		"There was an error in reading the data. Are you sure you've understood the data format?"
		<< std::endl;
	}
	else{
		std::cout << error_cause << std::endl;
	}
	std::exit(0);
}

void check_h5_read_call(herr_t ret){
	if (ret<0){
		raise_error("Erroneous control return-value from calling H5Dread()");
	}
}

void check_h5_write_call(herr_t ret){
	if (ret<0){
		raise_error("Erroneous control return-value from calling H5Dwrite()");
	}
}

void check_h5_read_table_call(herr_t ret){
	if (ret<0){
		raise_error("Erroneous control return-value from calling H5TBread_table()");
	}
}

void check_h5_close_call(herr_t ret){
	if (ret<0){
		raise_error("Erroneous control return-value from calling H5Dclose()");
	}
}
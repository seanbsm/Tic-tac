
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

bool verify_sparse_mkl_success(sparse_status_t mkl_sparse_status, bool print_true){
    bool operation_success = false;
    if (mkl_sparse_status==SPARSE_STATUS_SUCCESS){
        if (print_true){
            std::cout << "The operation was successful.                             " << std::endl;
        }
        operation_success = true;
    }
    else if (mkl_sparse_status==SPARSE_STATUS_NOT_INITIALIZED){
        std::cout << "The routine encountered an empty handle or matrix array.  " << std::endl;
    }
    else if (mkl_sparse_status==SPARSE_STATUS_ALLOC_FAILED){
        std::cout << "Internal memory allocation failed.                        " << std::endl;
    }
    else if (mkl_sparse_status==SPARSE_STATUS_INVALID_VALUE){
        std::cout << "The input parameters contain an invalid value.            " << std::endl;
    }
    else if (mkl_sparse_status==SPARSE_STATUS_EXECUTION_FAILED){
        std::cout << "Execution failed.                                         " << std::endl;
    }
    else if (mkl_sparse_status==SPARSE_STATUS_INTERNAL_ERROR){
        std::cout << "An error in algorithm implementation occurred.            " << std::endl;
    }
    else if (mkl_sparse_status==SPARSE_STATUS_NOT_SUPPORTED){
        std::cout << "The requested operation is not supported.                 " << std::endl;
    }
    else{
        std::cout << "Unknown value for mkl_sparse_status.                      " << std::endl;
    }
    return operation_success;
}
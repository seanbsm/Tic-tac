
#include "error_management.h"

void raise_error (std::string error_cause){
    if (error_cause == "problem_reading_data"){
        std::cout <<
        "There was an error in reading the data. Are you sure you've understood the data format?"
        << std::endl;
    }
}
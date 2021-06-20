
#ifndef SET_RUN_PARAMETERS_H
#define SET_RUN_PARAMETERS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cctype>
#include <algorithm>

#include "type_defs.h"
#include "store_functionality.h"

void read_input_list_and_set_parameters(run_params& run_parameters);
void show_usage();

void no_options_entered();

void use_input_list(run_params& run_parameters);

void set_default_values(run_params& run_parameters);

void unrecognised_option(std::string option);

void only_show_usage();

void set_run_parameters(int& argc, char* argv[], run_params& run_parameters);

void save_run_parameters_to_file(run_params& run_parameters);

#endif // SET_RUN_PARAMETERS_H

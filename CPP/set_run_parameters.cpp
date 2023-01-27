
#include "set_run_parameters.h"

/* Template soluation copied from: https://stackoverflow.com/questions/16605967/set-precision-of-stdto-string-when-converting-floating-point-values
 * on 21/02/2022.
 * to_string does not allow for adjusting number of printed decimals. This template is a workaround. */
template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6){
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

std::string type_to_string(bool input){
	std::string return_val = "false";
	if (input==true){
		return_val = "true";
	}
	return return_val;
}
std::string type_to_string(int input){
	return std::to_string(input);
}
std::string type_to_string(float input){
	return to_string_with_precision(input, 3);
}
std::string type_to_string(double input){
	return to_string_with_precision(input, 3);
}
std::string type_to_string(std::string input){
	return input;
}

std::string create_input_printout_string(run_params run_parameters){
	std::ostringstream output_string;
	output_string << "Running program for:" << std::endl;
	output_string << "two_J_3N_max:                    " << type_to_string(run_parameters.two_J_3N_max) 		  	  	<< "\n";
	output_string << "Np_WP:                           " << type_to_string(run_parameters.Np_WP) 				  	  	<< "\n";
	output_string << "Nq_WP:                           " << type_to_string(run_parameters.Nq_WP) 				  	  	<< "\n";
	output_string << "J_2N_max:                        " << type_to_string(run_parameters.J_2N_max) 			  	  	<< "\n";
	output_string << "Nphi:                            " << type_to_string(run_parameters.Nphi) 				  	  	<< "\n";
	output_string << "Nx:                              " << type_to_string(run_parameters.Nx) 				  	  		<< "\n";
	output_string << "chebyshev sparseness:            " << type_to_string(run_parameters.chebyshev_t) 		  	  		<< "\n";
	output_string << "chebyshev scale:                 " << type_to_string(run_parameters.chebyshev_s) 		  	  		<< "\n";
	output_string << "Np_per_WP:                       " << type_to_string(run_parameters.Np_per_WP) 			  	  	<< "\n";
	output_string << "Nq_per_WP:                       " << type_to_string(run_parameters.Nq_per_WP) 			  	  	<< "\n";
	output_string << "P123-recovery mode on:           " << type_to_string(run_parameters.P123_recovery) 		  	  	<< "\n";
	output_string << "P123 omp number of threads:      " << type_to_string(run_parameters.P123_omp_num_threads) 	  	<< "\n";
	output_string << "Tensor-force on:                 " << type_to_string(run_parameters.tensor_force) 		  	  	<< "\n";
	output_string << "Isospin-breaking in 1S0:         " << type_to_string(run_parameters.isospin_breaking_1S0) 	  	<< "\n";
	output_string << "Mid-point approximation:         " << type_to_string(run_parameters.midpoint_approx) 	  	  		<< "\n";
	output_string << "Calculate P123 and store:        " << type_to_string(run_parameters.calculate_and_store_P123) 	<< "\n";
	output_string << "Calculate breakup amplitudes:    " << type_to_string(run_parameters.include_breakup_channels) 	<< "\n";
	output_string << "Solve Faddeev equation:          " << type_to_string(run_parameters.solve_faddeev)   	  	  		<< "\n";
	output_string << "Solve Faddeev with LAPACK:       " << type_to_string(run_parameters.solve_dense)   	  	  		<< "\n";
	output_string << "Production run:                  " << type_to_string(run_parameters.production_run)  	  	  		<< "\n";
	output_string << "Potential model:                 " << type_to_string(run_parameters.potential_model) 	  	  		<< "\n";
	output_string << "p-momentum grid type:            " << type_to_string(run_parameters.p_grid_type) 		  	  		<< "\n";
	output_string << "p-momentum grid input file:      " << type_to_string(run_parameters.p_grid_filename) 	  	  		<< "\n";
	output_string << "q-momentum grid type:            " << type_to_string(run_parameters.q_grid_type) 		  	  		<< "\n";
	output_string << "q-momentum grid input file:      " << type_to_string(run_parameters.q_grid_filename) 	  	  		<< "\n";
	output_string << "Parameter walk:                  " << type_to_string(run_parameters.parameter_walk)		  		<< "\n";
	output_string << "Parameter input file:            " << type_to_string(run_parameters.parameter_file)		  		<< "\n";
	output_string << "Parameter set index (PSI) start: " << type_to_string(run_parameters.PSI_start)		  			<< "\n";
	output_string << "Parameter set index (PSI) end:   " << type_to_string(run_parameters.PSI_end)		  				<< "\n";
	output_string << "Output folder:                   " << type_to_string(run_parameters.output_folder) 		  		<< "\n";
	output_string << "P123-matrix read/write folder:   " << type_to_string(run_parameters.P123_folder) 		  	  		<< "\n";
	output_string << "Parallel run:                    " << type_to_string(run_parameters.parallel_run)  		  	  	<< "\n";
	if(run_parameters.parallel_run==true){
	output_string << "Channel index:                   " << type_to_string(run_parameters.channel_idx) 		  	  		<< "\n";
	}
	return output_string.str();
}

bool read_and_set_parameter(run_params& run_parameters, std::string option, std::string input){

	bool valid_option_found = true;
	
	if (option == "two_J_3N_max"){
		run_parameters.two_J_3N_max = std::stoi(input);
	}
	else if (option == "Np_WP"){
		run_parameters.Np_WP = std::stoi(input);
	}
	else if (option == "Nq_WP"){
		run_parameters.Nq_WP = std::stoi(input);
	}
	else if (option == "J_2N_max"){
		run_parameters.J_2N_max = std::stoi(input);
	}
	else if (option == "Nphi"){
		run_parameters.Nphi = std::stoi(input);
	}
	else if (option == "Nx"){
		run_parameters.Nx = std::stoi(input);
	}
	else if (option == "chebyshev_t"){
		run_parameters.chebyshev_t = std::stod(input);
	}
	else if (option == "chebyshev_s"){
		run_parameters.chebyshev_s = std::stod(input);
	}
	else if (option == "Np_per_WP"){
		run_parameters.Np_per_WP = std::stoi(input);
	}
	else if (option == "Nq_per_WP"){
		run_parameters.Nq_per_WP = std::stoi(input);
	}
	else if (option == "P123_omp_num_threads"){
		run_parameters.P123_omp_num_threads = std::stoi(input);
	}
	else if (option == "PSI_start"){
		run_parameters.PSI_start = std::stoi(input);
	}
	else if (option == "PSI_end"){
		run_parameters.PSI_end = std::stoi(input);
	}
	else if (option=="channel_idx"){
		run_parameters.channel_idx = std::stoi(input);
	}
	else if (option == "parallel_run"){
		if (input=="true" || input=="false"){
			run_parameters.parallel_run = (input=="true");
		}
		else{
			raise_error("Invalid value for input parameter parallel_run!");
		}
	}
	else if (option == "P123_recovery"){
		if (input=="true" || input=="false"){
			run_parameters.P123_recovery = (input=="true");
		}
		else{
			raise_error("Invalid value for input parameter P123_recovery!");
		}
	}
	else if (option == "tensor_force"){
		if (input=="true" || input=="false"){
			run_parameters.tensor_force = (input=="true");
		}
		else{
			raise_error("Invalid value for input parameter tensor_force!");
		}
	}
	else if (option == "isospin_breaking_1S0"){
		if (input=="true" || input=="false"){
			run_parameters.isospin_breaking_1S0 = (input=="true");
		}
		else{
			raise_error("Invalid value for input parameter isospin_breaking_1S0!");
		}
	}
	else if (option == "midpoint_approx"){
		if (input=="true" || input=="false"){
			run_parameters.midpoint_approx = (input=="true");
		}
		else{
			raise_error("Invalid value for input parameter midpoint_approx!");
		}
	}
	else if (option == "calculate_and_store_P123"){
		if (input=="true" || input=="false"){
			run_parameters.calculate_and_store_P123 = (input=="true");
		}
		else{
			raise_error("Invalid value for input parameter calculate_and_store_P123!");
		}
	}
	else if (option == "include_breakup_channels"){
		if (input=="true" || input=="false"){
			run_parameters.include_breakup_channels = (input=="true");
		}
		else{
			raise_error("Invalid value for input parameter include_breakup_channels!");
		}
	}
	else if (option == "solve_faddeev"){
		if (input=="true" || input=="false"){
			run_parameters.solve_faddeev = (input=="true");
		}
		else{
			raise_error("Invalid value for input parameter solve_faddeev!");
		}
	}
	else if (option == "solve_dense"){
		if (input=="true" || input=="false"){
			run_parameters.solve_dense = (input=="true");
		}
		else{
			raise_error("Invalid value for input parameter solve_faddeev!");
		}
	}
	else if (option == "production_run"){
		if (input=="true" || input=="false"){
			run_parameters.production_run = (input=="true");
		}
		else{
			raise_error("Invalid value for input parameter production_run!");
		}
	}
	else if (option == "parameter_walk"){
		if (input=="true" || input=="false"){
			run_parameters.parameter_walk = (input=="true");
		}
		else{
			raise_error("Invalid value for input parameter parameter_walk!");
		}
	}
	else if (option == "potential_model"){
		run_parameters.potential_model = input;
	}
	else if (option == "subfolder"){
		run_parameters.subfolder = input;
	}
	else if (option == "p_grid_type"){
		run_parameters.p_grid_type = input;
	}
	else if (option == "p_grid_filename"){
		run_parameters.p_grid_filename = input;
	}
	else if (option == "q_grid_type"){
		run_parameters.q_grid_type = input;
	}
	else if (option == "q_grid_filename"){
		run_parameters.q_grid_filename = input;
	}
	else if (option == "parameter_file"){
		run_parameters.parameter_file = input;
	}
	else if (option == "energy_input_file"){
		run_parameters.energy_input_file = input;
	}
	else if (option == "output_folder"){
		run_parameters.output_folder = input;
	}
	else if (option == "P123_folder"){
			run_parameters.P123_folder = input;
		}
	else{
		valid_option_found = false;
	}

	return valid_option_found;
}

void read_input_list_and_set_parameters(run_params& run_parameters, std::string filename){
	/* Define file with input */
	std::ifstream infile(filename);
	
	std::string line;
	
	std::string delimiter = "=";
	std::string option;
	std::string input;
	
	/* Loop through lines in file */
	while (std::getline(infile, line)){
		
		if (line.find(delimiter) != std::string::npos){
			
			option = line.substr(0, line.find(delimiter));
			input  = line.substr(line.find(delimiter)+1,-1);
			
			/* See if input is empty (only whitespace) */
			if (input.find_first_not_of(" ") == std::string::npos){
				continue;
			}
			
			/* Remove whitespace from input string */
			input.erase(std::remove_if(input.begin(), input.end(), ::isspace), input.end());
			
			bool valid_option_found = read_and_set_parameter(run_parameters, option, input);

			if (valid_option_found==false){
				unrecognised_option(line);
			}
		}
	}
}

void show_usage(){
	
	std::string seperationLine = "\n---------------------------------------------------------------------------------------\n";
	
	std::cout << seperationLine
			  << std::endl;

	std::cout << "NOTE: The code rewrites variables IN ORDER, meaning if a variable is submitted\n"
	          << "      twice (either on command line or in input-files), the last read variable\n"
			  << "      value will be used! \n"
			  << seperationLine
			  << std::endl;

	std::cout << "The following command-line run options are available:\n"
			  << seperationLine
			  << std::endl;
			  
	std::cout << "-h or --help:             Displays the possible run options.\n"
			  << seperationLine
			  << std::endl;

	std::cout << "<inputfile>.txt:          Tell program to read <inputfile>.txt for run parameters.\n"
			  << "Example:                  input.txt -> Program interprets input.txt.\n"
			  << seperationLine
			  << std::endl;

	std::cout << "two_J_3N_max:             Maximum total angular momentum of three-nucleon system state\n"
			  << "                          space \n"
			  << "Example:                  two_J_3N_max=3 -> J_3N=[1/2, 3/2].\n"
			  << seperationLine
			  << std::endl;
			  
	std::cout << "Np_WP:                    Sets number of p-momentum wave-packets.\n"
			  << "Example:                  N=100 -> 100 wave-packets in p-momentum.\n"
			  << seperationLine
			  << std::endl;
		
	std::cout << "Nq_WP:                    Sets number of q-momentum wave-packets.\n"
			  << "Example:                  N=100 -> 100 wave-packets in q-momentum.\n"
			  << seperationLine
			  << std::endl;

	std::cout << "J_2N_max:                 Sets upper bound on total angular momentum of\n"
			  << "                          pair-system. \n"
			  << "Example:                  J_2N_max=10 -> Calculates up to and including J_2N=10.\n"
			  << seperationLine
			  << std::endl;
	
	std::cout << "Nphi:                     Number of quadrature points in p- and q- parametrisation in\n"
			  << "                          P123-calculation.\n"
			  << "Example:                  Nphi=48 -> 48 quadrature points.\n"
			  << seperationLine
			  << std::endl;
	
	std::cout << "Nx:                       Number of quadrature points in angular integral of geometric\n"
			  << "                          function in P123-calculation.\n"
			  << "Example:                  Nx=20 -> 20 quadrature points.\n"
			  << seperationLine
			  << std::endl;
	
	std::cout << "chebyshev_t:              Sparseness degree in Chebyshev distribution.\n"
			  << "Example:                  t=2 -> Use exponent 2 in Chebyshev boundary construction.\n"
			  << seperationLine
			  << std::endl;
	
	std::cout << "chebyshev_s:              Scale in Chebyshev distribution.\n"
			  << "Example:                  s=100 -> Scale Chebyshev boundaries by 100 MeV.\n"
			  << seperationLine
			  << std::endl;
	
	std::cout << "Np_per_WP:                Sets number of quadrature points in each p-momentum wave\n"
			  << "                          packet. Used in NN potential matrix construction\n"
			  << "Example:                  Np_per_WP=8 -> 8 Gauss-Legendre points in each wave-packet\n"
			  << "                                         in p-momentum.\n"
			  << seperationLine
			  << std::endl;
		
	std::cout << "Nq_per_WP:                Sets number of quadrature points in each q-momentum wave\n"
	          << "                          packet.\n"
			  << "Example:                  Nq_per_WP=8 -> 8 Gauss-Legendre points in each wave-packet\n"
			  << "                                         in q-momentum.\n"
			  << seperationLine
			  << std::endl;

	std::cout << "channel_idx:              Tells the program to calculate for a given channel index\n"
			  << "                          (0-based indexing!).\n"
			  << "Example:                  channel_idx=2 -> Program calculates U/P123 for channel 2.\n"
			  << seperationLine
			  << std::endl;

	std::cout << "parallel_run:             Tells the program to run in parallel across 3N channels (one\n"
			  << "                          processor per channel). Possible options are (true, false).\n"
			  << "Example:                  parallel_run=true -> Program assigns channels to processors\n"
			  << "                                               given by channel_idx.\n"
			  << seperationLine
			  << std::endl;

	std::cout << "potential_model:          Sets which two-body potential model is used. Possible\n"
			  << "                          options are (LO_internal, N2LOopt, Idaho_N3LO, nijmegen,\n"
			  << "                          malfliet_tjon). Note that LO_internal is an internally pre-\n"
			  << "                          written chiral leading-order potential.\n"
			  << "Example:                  model=Idaho_N3LO -> Uses two-body chiral Idaho N3LO\n"
	  		  << "                                              potential.\n"
			  << seperationLine
			  << std::endl;
	
	std::cout << "p_grid_type:              Sets which type of distribution of p-momentum bins.\n"
			  << "                          Currently only the Chebyshev distribution is implemented.\n"
			  << "Example:                  p_grid_type=chebyshev -> WP-Boundaries are distributed in a\n"
			  << "                                                   Chebyshev-distribution.\n"
			  << seperationLine
			  << std::endl;
	
	std::cout << "p_grid_filename:          Input-file with custom-set boundaries of p-momentum bins.\n"
			  << "Example:                  p_grid_filename=filename.txt -> Program reads boundaries in\n"
			  << "                                                          filename.txt. Must be stored\n"
			  << "                                                          in rows.\n"
			  << seperationLine
			  << std::endl;
	
	std::cout << "q_grid_type:              Same as p_grid_type but for q-momentum bins.\n"
			  << seperationLine
			  << std::endl;
	
	std::cout << "q_grid_filename:          Same as p_grid_filename but for q-momentum bins.\n"
			  << seperationLine
			  << std::endl;
			  
	std::cout << "P123_omp_num_threads:     Number of OpenMP threads to use in calculation of P123-\n"
	          << "                          elements\n"
			  << "Example:                  P123_omp_num_threads=8 -> Calculation uses 8 threads.\n"
			  << seperationLine
			  << std::endl;
			  
	std::cout << "parameter_walk:           Sets whether or not program is run as a walk over model\n"
			  << "                          parameter space. (REQUIRES INPUT FILE, SEE parameter_file)\n"
			  << "                          Possible options are (true, false).\n"
			  << "Example:                  parameter_walk=true -> Program loops over one or several \n"
			  << "                                                 model parameter sets.\n"
			  << seperationLine
			  << std::endl;
			  
	std::cout << "parameter_file:           File from which the code reads parameter input. Note that\n"
			  << "                          the number of parameters must equal the number of potential\n"
	  		  << "                          model parameters \n"
			  << "Example:                  parameter_walk=Input/samples.txt -> Program reads samples.txt\n"
	  		  << "                                                              line-by-line and loops\n"
			  << "                                                              through values. \n"
			  << seperationLine
			  << std::endl;
			  
	std::cout << "PSI_start:                Index of parameter set to start iteration through\n"
			  << "                          parameter_file. (!!! ZERO-BASED INDEXING !!!) \n"
			  << "Example:                  PSI_start=10 -> Program reads parameter_file, and starts\n"
	  		  << "                          iteration of parameters sets at line 10. \n"
			  << seperationLine
			  << std::endl;
			  
	std::cout << "PSI_end:                  Index of parameter set to end iteration through\n"
			  << "                          parameter_file. (!!! ZERO-BASED INDEXING !!!) \n"
			  << "Example:                  PSI_end=20 -> Program reads parameter_file, and ends\n"
	  		  << "                                        iteration of parameters sets at line 19 (not\n"
	  		  << "                                        using line 20). Combined example with\n"
	  		  << "                                        PSI_start: if PSI_start=0 and PSI_end=10, the\n"
	  		  << "                                        code requires a file with at least 10 lines of\n"
			  << "                                        parameter sets, where the first line is given\n"
			  << "                                        index 0, such that PSI=[0,1,...,9]."
			  << seperationLine
			  << std::endl;
			  
	std::cout << "P123_recovery:            \n"
			  << seperationLine
			  << std::endl;
			  
	std::cout << "tensor_force:             Tell program whether to include tensor forces in nuclear NN\n"
			  << "                          potential. Possible options are (true, false). \n"
			  << "Example:                  tensor_force=true -> Tensor forces are enabled.\n"
			  << seperationLine
			  << std::endl;
			  
	std::cout << "isospin_breaking_1S0:     Tell program whether to include charge depence in the 1S0 NN\n"
			  << "                          channel. Possible options are (true, false). \n"
			  << "Example:                  isospin_breaking_1S0=true -> charge dependence is included.\n"
			  << seperationLine
			  << std::endl;
	
	std::cout << "midpoint_approx:          Sets whether or not program will use momentum bin averages\n"
			  << "                          when calculating potential elements. Possible options are\n"
	  		  << "                          (true, false).\n"
			  << "Example:                  midpoint_approx=true -> Program will use average momentum of\n"
	  		  << "                                                  bins.\n"
			  << seperationLine
			  << std::endl;
			  
	std::cout << "calculate_and_store_P123: Tell program whether to calculate P123-elements and store.\n"
			  << "                          Possible options are (true, false).\n"
			  << "Example:                  calculate_and_store_P123=true -> Program will calculate P123\n"
	  		  << "                                                           elements and store them.\n"
			  << seperationLine
			  << std::endl;
			  
	std::cout << "include_breakup_channels: Tell program whether to calculate breakup channels, This is\n"
			  << "                          more costly. Possible options are (true, false).\n"
			  << "Example:                  include_breakup_channels=true -> Program will also calculate \n"
	  		  << "                                                           breakup U-matrix elements.\n"
			  << seperationLine
			  << std::endl;
			  
	std::cout << "solve_faddeev:            Tell program whether to solve the Faddeev equation.\n"
			  << "                          Possible options are (true, false).\n"
			  << "Example:                  solve_faddeev=true -> Program will solve the Faddeev\n"
	  		  << "                                                equation.\n"
			  << seperationLine
			  << std::endl;
			  
	std::cout << "solve_dense:              Tell program whether to solve the Faddeev equation using\n"
			  << "                          a dense LAPACK solver or with Pade resummation. Note that\n"
			  << "                          the dense solver is VERY SLOW compared to the Pade resummation.\n"
			  << "                          The dense solver also stores the A, K, and U matrices. This\n"
			  << "                          option only really exists for debugging.\n"
			  << "                          Possible options are (true, false).\n"
			  << "Example:                  solve_dense=true -> Program will solve the Faddeev equation\n"
	  		  << "                                              using LAPACK solver.\n"
			  << seperationLine
			  << std::endl;
			  
	std::cout << "production_run:           Tell program to do correct calculations or not. Used to\n"
			  << "                          debug P123-calculations quickly with nonsense elements\n"
			  << "                          Possible options are (true, false).\n"
			  << "Example:                  production_run=false -> Program will calculate P123\n"
	  		  << "                                                  with low-cost placeholder elements\n"
			  << "                                                  in sparse layout.\n"
			  << seperationLine
			  << std::endl;
	
	std::cout << "energy_input_file:        Sets in which file in Input to read for Tlab energies used\n"
			  << "                          in U-matrix calculation.\n"
			  << "Example:                  energy_input_file=Tlabs -> Program reads file <Tlabs.txt> \n"
			  << seperationLine
			  << std::endl;
	
	std::cout << "output_folder:            Sets in which folder store output into. \n"
			  << "Example:                  output_folder=some/folder -> Program stores results in \n"
			  << "                                                       some/folder/ \n"
			  << seperationLine
			  << std::endl;
	
	std::cout << "P123_folder:              Sets in which folder to read/write P123 from/to.\n"
			  << "Example:                  P123_folder=some/folder -> Program reads matrix from, or\n"
			  << "                                                     writes matrix to some/folder/ \n"
			  << seperationLine
			  << std::endl;
	
	std::cout << std::endl;
}

void no_options_entered(){
	std::cout << "No run options entered." << std::endl;
	std::cout << "Use -h or --help for run options." << std::endl;
	std::cout << std::endl;
	std::cout << "Running default setup." << std::endl;
	std::cout << std::endl;
}

void use_input_list(run_params& run_parameters, std::string filename){
						
	read_input_list_and_set_parameters(run_parameters, filename);
	
	std::cout << "Using input file " <<  filename << " for run options." << std::endl;
	std::cout << std::endl;
	std::cout << "Any unset options will be filled with default parameters." << std::endl;
	std::cout << "All other command line input will be ignored." << std::endl;
	std::cout << std::endl;
}

void unrecognised_option(std::string option){
	std::cout << "Unrecognised run option entered: " << option << std::endl;
	std::cout << "Use -h or --help for run options." << std::endl;
	std::cout << std::endl;
	std::cout << "Skipping unrecognised run option." << std::endl;
	std::cout << std::endl;
}

void only_show_usage(){
	std::cout << "Only -h or --help was requested. " << std::endl;
	std::cout << "Program will not run. Exiting ..." << std::endl;
	std::cout << std::endl;
}

void set_default_values(run_params& run_parameters){

	/* Default parameters */
	run_parameters.two_J_3N_max 	        = 1;
	run_parameters.Np_WP		 	        = 50;
	run_parameters.Nq_WP		 	        = 50;
	run_parameters.J_2N_max	  	  	        = 3;
	run_parameters.Nphi		 		        = 48;
	run_parameters.Nx 			 	        = 48;
	run_parameters.chebyshev_t		        = 1;
	run_parameters.chebyshev_s		        = 100;
	run_parameters.Np_per_WP	 	        = 8;
	run_parameters.Nq_per_WP	 	        = 8;
	run_parameters.channel_idx		        = -1;
	run_parameters.parallel_run		        = false;
	run_parameters.potential_model	        = "LO_internal";
	run_parameters.subfolder	  	        = "Output";
	run_parameters.p_grid_type 	  	        = "chebyshev";
	run_parameters.p_grid_filename 	        = "";
	run_parameters.q_grid_type 	  	        = "chebyshev";
	run_parameters.q_grid_filename 	        = "";
	run_parameters.P123_omp_num_threads     = omp_get_max_threads();
	run_parameters.parameter_walk 	        = false;
	run_parameters.parameter_file 	        = "none";
	run_parameters.PSI_start				= -1;
	run_parameters.PSI_end				    = -1;
	run_parameters.P123_recovery		    = false;
	run_parameters.tensor_force			    = true;
	run_parameters.isospin_breaking_1S0     = true;
	run_parameters.midpoint_approx  		= false;
	run_parameters.calculate_and_store_P123 = true;
	run_parameters.include_breakup_channels = false;
	run_parameters.solve_faddeev		    = true;
	run_parameters.solve_dense				= false;
	run_parameters.production_run		    = true;
	run_parameters.energy_input_file        = "lab_energies.txt";
	run_parameters.output_folder  	        = "Output";
	run_parameters.P123_folder  	        = "Output";
}

void set_run_parameters(int& argc, char* argv[], run_params& run_parameters){
	
	set_default_values(run_parameters);
	
	if (argc == 1) {
		no_options_entered();
	}
	else{
		std::string delimiter = "=";
		std::string arg;
		std::string option;
		std::string input;
		
		for (int i = 1; i < argc; ++i) {
			
			arg = argv[i];
			size_t arg_size = arg.size();
			
			if ((arg == "-h") || (arg == "--help")) {
				show_usage();
				if (argc == 2) {
					only_show_usage();
					exit(-1);
				}
			}
			else if (arg.substr(arg_size-4) == ".txt"){
				/* Read input list and set input variables specified */
				use_input_list(run_parameters, arg);
			}
			else if (arg.find(delimiter) != std::string::npos){
				
				option = arg.substr(0, arg.find(delimiter));
				input  = arg.substr(arg.find(delimiter)+1,-1);
				
				if (input==""){
					continue;
				}
				
				bool valid_option_found = read_and_set_parameter(run_parameters, option, input);

				if (valid_option_found==false){
					unrecognised_option(arg);
				}
			}
			else{
				unrecognised_option(arg);
			}
		}
	}

	if (run_parameters.parallel_run==true && run_parameters.channel_idx==-1){
		raise_error("Parallel run specified in input but no channel index (i.e. channel_idx=<int>) given!");
	}

	/* Do program compatibility-checks with input */
	if (run_parameters.P123_omp_num_threads>run_parameters.Nq_WP){
		run_parameters.P123_omp_num_threads = run_parameters.Nq_WP;
		std::cout << "NOTE: Nq_WP is smaller than number of threads. P123 parallellism is with respect to Nq_WP." << std::endl;
		std::cout << "      -> Setting P123_omp_num_threads=Nq_WP." << std::endl;
	}
	if ( run_parameters.two_J_3N_max%2==0 ||  run_parameters.two_J_3N_max<=0 ){
		raise_error("Cannot have even two_J_3N_max!");
	}

	/* Print system run parameters */
	std::cout << std::endl;
	std::cout << create_input_printout_string(run_parameters);
	std::cout << std::endl;
	
	store_run_parameters(run_parameters);
}

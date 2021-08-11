
#include "set_run_parameters.h"

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
			else if (option == "parallel_run"){
				if (input=="true" || input=="false"){
					run_parameters.parallel_run = (input=="true");
				}
				else{
					raise_error("Invalid value for input parameter parallel_run!");
				}
			}
			else if (option == "potential_model"){
				run_parameters.potential_model = input;
			}
			else if (option == "subfolder"){
				run_parameters.subfolder = input;
			}
			else if (option == "grid_type"){
				run_parameters.grid_type = input;
			}
			else if (option == "parameter_walk"){
				run_parameters.parameter_walk = input;
			}
			else if (option == "energy_input_file"){
				run_parameters.energy_input_file = input;
			}
			else if (option == "average"){
				run_parameters.average = input;
			}
			else if (option == "output_folder"){
				run_parameters.output_folder = input;
			}
			else if (option == "P123_folder"){
					run_parameters.P123_folder = input;
				}
			else{
				unrecognised_option(line);
			}
		}
	}
}

void show_usage(){
	
	std::string seperationLine = "\n --------------------------------------------------------------------------- \n";
	
	std::cout << "The following command-line run options are available:\n"
			  << seperationLine
			  << std::endl;
			  
	std::cout << "-h or --help: Displays the possible run options.\n"
			  << seperationLine
			  << std::endl;

	std::cout << "<inputfile>.txt:    Tell program to read <inputfile>.txt for run parameters. \n"
			  << "Example:            input.txt -> Program interprets input.txt.\n"
			  << "Example:                         All other command-line options are ignored.\n"
			  << seperationLine
			  << std::endl;

	std::cout << "two_J_3N_max:       Maximum total angular momentum of three-nucleon system state space \n"
			  << "Example:            two_J_3N_max=3 -> J_3N=[1/2, 3/2].\n"
			  << seperationLine
			  << std::endl;
			  
	std::cout << "Np_WP:              Sets number of p-momentum wave-packets.\n"
			  << "Example:            N=100 -> 100 wave-packets in p-momentum.\n"
			  << seperationLine
			  << std::endl;
		
	std::cout << "Nq_WP:              Sets number of q-momentum wave-packets.\n"
			  << "Example:            N=100 -> 100 wave-packets in q-momentum.\n"
			  << seperationLine
			  << std::endl;

	std::cout << "J_2N_max:           Sets upper bound on total angular momentum of\n"
			  << "                    pair-system. \n"
			  << "Example:            J_2N_max=10 -> Calculates up to and including J_2N=10.\n"
			  << seperationLine
			  << std::endl;
	
	std::cout << "Nphi:               Number of quadrature points in p- and q- parametrisation in P123-calculation.\n"
			  << "Example:            Nphi=48 -> 48 quadrature points.\n"
			  << seperationLine
			  << std::endl;
	
	std::cout << "Nx:                 Number of quadrature points in angular integral of geometric function in P123-calculation.\n"
			  << "Example:            Nx=20 -> 20 quadrature points.\n"
			  << seperationLine
			  << std::endl;

	std::cout << "channel_idx:        Tells the program to calculate for a given channel index (0-based indexing!).\n"
			  << "Example:            channel_idx=2 -> Program calculates U/P123 for channel 2.\n"
			  << seperationLine
			  << std::endl;

	std::cout << "parallel_run:       Tells the program to run in parallel across 3N channels (one processor per channel).\n"
			  << "                    Possible options are (true, false).\n"
			  << "Example:            parallel_run=true -> Program assigns channels to processors given by channel_idx.\n"
			  << seperationLine
			  << std::endl;

	std::cout << "potential_model:    Sets which two-bodu potential model is used.\n"
			  << "                    Possible options are (LO_internal, N2LOopt, Idaho_N3LO, nijmegen, malfliet_tjon).\n"
			  << "                    Note that LO_internal is an internally pre-written chiral leading-order potential.\n"
			  << "Example:            model=Idaho_N3LO -> Uses two-body chiral Idaho N3LO potential.\n"
			  << seperationLine
			  << std::endl;
	
	std::cout << "grid_type:          Sets which type of distribution of bins.\n"
			  << "                    Currently only the Chebyshev distribution is implemented.\n"
			  << "Example:            grid_type=chebyshev -> WP-Boundaries are distributed in a Chebyshev-distribution.\n"
			  << seperationLine
			  << std::endl;
			  
	std::cout << "parameter_walk:     Sets whether or not program is run as a walk over \n"
			  << "                    model parameter space.\n"
			  << "                    Possible options are (on, off).\n"
			  << "Example:            parameter_walk=on -> Program loops over several \n"
			  << "				                           model parameter sets.\n"
			  << seperationLine
			  << std::endl;
	
	std::cout << "energy_input_file:  Sets in which file in Input to read for Tlab energies used in T-matrix calculation.\n"
			  << "Example:            energy_input_file=Tlabs -> Program reads file <Tlabs.txt> \n"
			  << seperationLine
			  << std::endl;
	
	std::cout << "average:            Sets whether or not program will use momentum bin averages \n"
			  << "                    when calculating potential elements.\n"
			  << "                    Possible options are (on, off).\n"
			  << "Example:            average=on -> Program will use average momentum of bins \n"
			  << seperationLine
			  << std::endl;
	
	std::cout << "output_folder:      Sets in which folder store output into. \n"
			  << "Example:            output_folder=some/folder -> Program stores results in \n"
			  << "				                                   some/folder/ \n"
			  << seperationLine
			  << std::endl;
	
	std::cout << "P123_folder:        Sets in which folder to read/write P123 from/to. \n"
			  << "Example:            P123_folder=some/folder -> Program reads matrix from, or writes matrix to, \n"
			  << "				                                 some/folder/ \n"
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
	run_parameters.two_J_3N_max 	  = 1;
	run_parameters.Np_WP		 	  = 50;
	run_parameters.Nq_WP		 	  = 50;
	run_parameters.J_2N_max	  	  	  = 3;
	run_parameters.Nphi		 		  = 48;
	run_parameters.Nx 			 	  = 15;
	run_parameters.Np_per_WP	 	  = 8;
	run_parameters.Nq_per_WP	 	  = 8;
	run_parameters.channel_idx		  = -1;
	run_parameters.parallel_run		  = false;
	run_parameters.potential_model	  = "LO_internal";
	run_parameters.subfolder	  	  = "Output";
	run_parameters.grid_type 	  	  = "chebyshev";
	run_parameters.parameter_walk 	  = "off";
	run_parameters.energy_input_file  = "lab_energies";
	run_parameters.average  	  	  = "off";
	run_parameters.output_folder  	  = "Output";
	run_parameters.P123_folder  	  = "Output";
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
			else if (arg.substr(arg_size-4) == ".txt" && run_parameters.parallel_run==false){
				/* Reset to default values.
				 * Essentially we ignore all other
				 * command line input. */
				set_default_values(run_parameters);
				
				use_input_list(run_parameters, arg);
				if (run_parameters.parallel_run==false){
					break;
				}
				else{
					for (int i = 1; i < argc; ++i) {
						arg = argv[i];
						size_t arg_size = arg.size();
						if (arg.find(delimiter) != std::string::npos){
				
							option = arg.substr(0, arg.find(delimiter));
							input  = arg.substr(arg.find(delimiter)+1,-1);

							if (input==""){
								continue;
							}
							else if (option=="channel_idx"){
								run_parameters.channel_idx = std::stoi(input);
							}
						}
					}
					if (run_parameters.channel_idx==-1){
						raise_error("Parallel run specified in input but no channel index given!");
					}
					else{
						break;
					}
				}
			}
			else if (arg.find(delimiter) != std::string::npos){
				
				option = arg.substr(0, arg.find(delimiter));
				input  = arg.substr(arg.find(delimiter)+1,-1);
				
				if (input==""){
					continue;
				}
				
				
				if (option == "two_J_3N_max"){
					run_parameters.two_J_3N_max= std::stoi(input);
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
				else if (option == "parallel_run"){
					if (input=="true" || input=="false"){
						run_parameters.parallel_run = (input=="true");
					}
					else{
						raise_error("Invalid value for input parameter parallel_run!");
					}
				}
				else if (option == "potential_model"){
					run_parameters.potential_model = input;
				}
				else if (option == "grid_type"){
					run_parameters.grid_type = input;
				}
				else if (option == "parameter_walk"){
					run_parameters.parameter_walk = input;
				}
				else if (option == "energy_input_file"){
					run_parameters.energy_input_file = input;
				}
				else if (option == "average"){
					run_parameters.average = input;
				}
				else if (option == "output_folder"){
					run_parameters.output_folder = input;
				}
				else if (option == "P123_folder"){
					run_parameters.P123_folder = input;
				}
				else{
					unrecognised_option(arg);
				}
			}
			else{
				unrecognised_option(arg);
			}
		}
	}

	/* Print system run parameters */
	std::cout << std::endl;
	std::cout << "Running program for:" << std::endl;
	std::cout << "two_J_3N_max:                  " << run_parameters.two_J_3N_max    << std::endl;
	std::cout << "Np_WP:                         " << run_parameters.Np_WP           << std::endl;
	std::cout << "Nq_WP:                         " << run_parameters.Nq_WP           << std::endl;
	std::cout << "J_2N_max:                      " << run_parameters.J_2N_max        << std::endl;
	std::cout << "Nphi:                          " << run_parameters.Nphi            << std::endl;
	std::cout << "Nx:                            " << run_parameters.Nx              << std::endl;
	if(run_parameters.parallel_run==true){
	std::cout << "Channel index:                 " << run_parameters.channel_idx     << std::endl;
	}
	std::cout << "Parallel run:                  " << run_parameters.parallel_run    << std::endl;
	std::cout << "Potential model:               " << run_parameters.potential_model << std::endl;
	std::cout << "Grid type:                     " << run_parameters.grid_type       << std::endl;
	std::cout << "Parameter walk:                " << run_parameters.parameter_walk  << std::endl;
	std::cout << "Bin averaging:                 " << run_parameters.average         << std::endl;
	std::cout << "Output folder:                 " << run_parameters.output_folder   << std::endl;
	std::cout << "P123-matrix read/write folder: " << run_parameters.P123_folder     << std::endl;
	std::cout << std::endl;
	
	store_run_parameters(run_parameters);
}
